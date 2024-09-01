import torch
from torch import nn
from transformers import (
    AutoTokenizer,
    AutoModelForSequenceClassification,
    Trainer,
    TrainingArguments,
)
from datasets import Dataset
from peft import get_peft_model, LoraConfig, TaskType
from sklearn.metrics import (
    precision_score,
    recall_score,
    f1_score,
    roc_auc_score,
    average_precision_score,
    confusion_matrix,
    matthews_corrcoef,
)
import numpy as np
import pandas as pd

MAX_SEQ_LENGTH = 151

tokeniser = AutoTokenizer.from_pretrained(
    "InstaDeepAI/nucleotide-transformer-2.5b-1000g"
)


class MethylationAwareModel(nn.Module):
    def __init__(self, base_model, seq_len=151):  # Ensure seq_len is correctly set
        super().__init__()
        self.base_model = base_model
        self.config = base_model.config
        self.seq_len = seq_len
        # Methylation Embedding for one-hot encoded input
        self.methylation_embedding = nn.Linear(4, base_model.config.hidden_size)
        # Combine Methylation and Sequence embeddings
        combined_hidden_size = base_model.config.hidden_size * 2
        # Classifier Head
        self.classifier = nn.Sequential(
            nn.Linear(combined_hidden_size, 512),
            nn.ReLU(),
            nn.Dropout(0.1),
            nn.Linear(512, 2),  # Binary classification
        )

    def forward(
        self,
        input_ids=None,
        attention_mask=None,
        methylation=None,
        labels=None,
        **kwargs,
    ):
        print(f"Model is in {'training' if self.training else 'evaluation'} mode.")
        print(f"Shape of input_ids: {input_ids.shape}")

        # Base model forward pass
        outputs = self.base_model.esm(
            input_ids=input_ids,
            attention_mask=attention_mask,
            output_hidden_states=True,
            return_dict=True,
        )
        sequence_output = outputs.last_hidden_state

        # Log the shape of sequence_output
        print(
            f"Shape of sequence_output during {'training' if self.training else 'evaluation'}: {sequence_output.shape}"
        )

        # Adjust methylation values
        methylation = methylation + 1  # Adjust to non-negative values
        methylation_one_hot = nn.functional.one_hot(
            methylation.long(), num_classes=4
        ).float()
        # Methylation embedding
        methylation_embeds = self.methylation_embedding(methylation_one_hot)
        # Log the shape of methylation_embeds
        print(
            f"Shape of methylation_embeds before concatenation: {methylation_embeds.shape}"
        )
        # Ensure consistency in sequence length
        assert sequence_output.size(1) == methylation_embeds.size(
            1
        ), f"Mismatch in sequence lengths: {sequence_output.size(1)} vs {methylation_embeds.size(1)}"
        # Concatenate sequence and methylation embeddings
        combined_output = torch.cat((sequence_output, methylation_embeds), dim=-1)
        # Log the shape of combined_output
        print(f"Shape of combined_output: {combined_output.shape}")
        logits = self.classifier(
            combined_output
        )  # Pass combined embeddings through classifier
        masked_logits = logits * attention_mask.unsqueeze(-1)
        pooled_logits = masked_logits.sum(dim=1) / attention_mask.sum(
            dim=1, keepdim=True
        )
        # Compute loss if labels are provided
        loss = None
        if labels is not None:
            loss_fct = nn.CrossEntropyLoss()
            loss = loss_fct(
                pooled_logits.view(-1, self.config.num_labels), labels.view(-1)
            )
        return (
            {"loss": loss, "logits": pooled_logits}
            if loss is not None
            else {"logits": pooled_logits}
        )


def process_methylation(sequence, base_methylation):
    threemer_methylation = []
    for i in range(1, len(sequence) - 1):
        if sequence[i] != "C":
            threemer_methylation.append(2)
        else:
            threemer_methylation.append(base_methylation[i])
    return threemer_methylation

def tokenize_function(examples):
    print(f"Original sequence: {examples['sequence'][0]}")  # Print the first sequence
    print(f"Length of original sequence: {len(examples['sequence'][0])}")
    
    outputs = tokeniser(
        examples["sequence"],
        truncation=False,  # Ensure no truncation
        padding=False,  # No padding; sequences should already be the correct length
        max_length=151,  # Ensure max_length is set to 151
    )
    
    print(f"Tokenized input_ids: {outputs['input_ids'][0]}")
    print(f"Length of tokenized input_ids: {len(outputs['input_ids'][0])}")
    
    processed_methylation = []
    for seq, methyl in zip(examples["sequence"], examples["methylation"]):
        processed_methylation.append(methyl)
    
    outputs["methylation"] = processed_methylation
    return outputs

def prepare_dataset(sequences, labels, methylation_data):
    dataset = Dataset.from_dict(
        {"sequence": sequences, "labels": labels, "methylation": methylation_data}
    )
    return dataset.map(tokenize_function, batched=True, remove_columns=["sequence"])


def compute_metrics(eval_pred):
    logits, labels = eval_pred
    predictions = np.argmax(logits, axis=-1)
    precision = precision_score(labels, predictions, average="binary")
    recall = recall_score(labels, predictions, average="binary")
    f1 = f1_score(labels, predictions, average="binary")
    accuracy = np.mean(predictions == labels)
    # AUC requires probability scores, which are logits here
    auc_roc = roc_auc_score(
        labels, logits[:, 1]
    )  # Assuming positive class is at index 1
    auc_pr = average_precision_score(labels, logits[:, 1])
    # MCC
    mcc = matthews_corrcoef(labels, predictions)
    # Confusion Matrix
    cm = confusion_matrix(labels, predictions)
    d = {
        "accuracy": accuracy,
        "precision": precision,
        "recall": recall,
        "f1_score": f1,
        "auc_roc": auc_roc,
        "auc_pr": auc_pr,
        "mcc": mcc,
        "confusion_matrix": cm.tolist(),  # Convert to list to avoid serialization issues
    }
    print(d)
    return d


def main():
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Using device: {device}")
    # Load datasets
    train_data = pd.read_parquet(
        "/users/zetzioni/sharedscratch/tapsformer/data/samples/methylation/raw_with_liver/top_250/train_data.parquet"
    )
    validation_data = pd.read_parquet(
        "/users/zetzioni/sharedscratch/tapsformer/data/samples/methylation/raw_with_liver/top_250/validation_data.parquet"
    )
    test_data = pd.read_parquet(
        "/users/zetzioni/sharedscratch/tapsformer/data/samples/methylation/raw_with_liver/top_250/test_data.parquet"
    )
    train_data["len"] = train_data.reference.str.len()
    validation_data["len"] = validation_data.reference.str.len()
    test_data["len"] = test_data.reference.str.len()
    train_data = train_data[train_data["len"] == 151]
    validation_data = validation_data[validation_data["len"] == 151]
    test_data = test_data[test_data["len"] == 151]
    train_sequences, train_labels, train_methylation = (
        train_data.reference.values,
        train_data.label.values,
        train_data.methylation_state.values,
    )
    val_sequences, val_labels, val_methylation = (
        validation_data.reference.values,
        validation_data.label.values,
        validation_data.methylation_state.values,
    )
    test_sequences, test_labels, test_methylation = (
        test_data.reference.values,
        test_data.label.values,
        test_data.methylation_state.values,
    )

    # Setup model
    base_model = AutoModelForSequenceClassification.from_pretrained(
        "InstaDeepAI/nucleotide-transformer-2.5b-1000g", num_labels=2
    )
    model = MethylationAwareModel(base_model)
    model = model.to(device)
    peft_config = LoraConfig(
        task_type=TaskType.SEQ_CLS,
        inference_mode=False,
        r=1,
        lora_alpha=32,
        lora_dropout=0.1,
        target_modules=["query", "value"],
    )
    lora_model = get_peft_model(model, peft_config)
    lora_model.print_trainable_parameters()
    lora_model = lora_model.to(device)
    train_dataset = prepare_dataset(train_sequences, train_labels, train_methylation)
    val_dataset = prepare_dataset(val_sequences, val_labels, val_methylation)
    test_dataset = prepare_dataset(test_sequences, test_labels, test_methylation)

    batch_size = 8
    training_args = TrainingArguments(
        f"tapsformer-finetuned-lora-NucleotideTransformer",
        remove_unused_columns=False,
        evaluation_strategy="steps",
        save_strategy="steps",
        learning_rate=5e-4,
        per_device_train_batch_size=batch_size,
        gradient_accumulation_steps=1,
        per_device_eval_batch_size=64,
        num_train_epochs=2,
        warmup_steps=500,
        weight_decay=0.01,
        logging_dir="./logs",
        logging_steps=100,
        eval_steps=100,
        save_steps=100,
        load_best_model_at_end=True,
        metric_for_best_model="f1_score",
        label_names=["labels"],
        dataloader_drop_last=True,
        max_steps=1000,
    )
    trainer = Trainer(
        # model=lora_model,
        model=lora_model,
        args=training_args,
        train_dataset=train_dataset,
        eval_dataset=val_dataset,
        tokenizer=tokeniser,
        compute_metrics=compute_metrics,
    )
    # Start training
    trainer.train()

    # Evaluate on test set
    test_results = trainer.evaluate(test_dataset)
    print(f"Test results: {test_results}")

    # Save the model
    trainer.save_model("./fine_tuned_methylation_model")


if __name__ == "__main__":
    main()
