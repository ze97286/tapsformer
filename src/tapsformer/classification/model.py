import torch.nn as nn
import torch.nn.functional as F
import torch

from transformers import AutoModel, AutoTokenizer

import numpy as np
import math, os
from copy import deepcopy

class ReadClassifier(nn.Module):
    def __init__(self, pretrained_model_name="InstaDeepAI/nucleotide-transformer-500m-human-ref", num_dmrs=100, seq_len=512):
        super().__init__()
        self.bert = AutoModel.from_pretrained(pretrained_model_name)
        self.tokenizer = AutoTokenizer.from_pretrained(pretrained_model_name)
        
        self.seq_len = seq_len
        
        # DMR embedding
        self.dmr_encoder = nn.Embedding(num_dmrs, self.bert.config.hidden_size)
        
        # Combine all features
        self.feature_combiner = nn.Linear(self.bert.config.hidden_size * 2, self.bert.config.hidden_size)
        
        # Improved classifier head
        self.classifier = nn.Sequential(
            nn.Linear(self.bert.config.hidden_size, 256),
            nn.ReLU(),
            nn.Dropout(0.1),
            nn.LayerNorm(256),
            nn.Linear(256, 64),
            nn.ReLU(),
            nn.Dropout(0.1),
            nn.LayerNorm(64),
            nn.Linear(64, 2)  # 2 classes: tumor vs healthy
        )

    def forward(self, input_ids, attention_mask, methylation_status, position, dmr_id):
        # Process DNA sequence through BERT
        bert_output = self.bert(input_ids=input_ids, attention_mask=attention_mask)
        sequence_output = bert_output.last_hidden_state
        
        # Encode DMR and add to sequence output
        dmr_encoded = self.dmr_encoder(dmr_id).unsqueeze(1).expand(-1, sequence_output.size(1), -1)
        combined_output = torch.cat([sequence_output, dmr_encoded], dim=-1)
        
        # Combine features
        combined_features = self.feature_combiner(combined_output)
        
        # Pool the sequence (mean pooling)
        pooled_output = combined_features.mean(dim=1)
        
        # Classify
        logits = self.classifier(pooled_output)
        
        return {
            "logits": logits,
            "dmr_logits": combined_features,
            "classification_logits": torch.softmax(logits, dim=-1)
        }

# Example usage
def prepare_input(dna_sequence, methylation_status, position, dmr_id, tokenizer):
    inputs = tokenizer(dna_sequence, return_tensors="pt", max_length=512, truncation=True, padding="max_length")
    methylation_tensor = torch.tensor(methylation_status, dtype=torch.float32)
    position_tensor = torch.tensor(position, dtype=torch.float32)
    dmr_tensor = torch.tensor(dmr_id, dtype=torch.long)
    return inputs["input_ids"], inputs["attention_mask"], methylation_tensor, position_tensor, dmr_tensor

# Instantiate the model
model = Classifier()

# Example input (you would need to prepare your actual input data)
dna_sequence = "ATCG" * 128  # Example sequence
methylation_status = [0, 1, 0, 1] * 128  # Example methylation status
position = list(range(512))  # Example position
dmr_id = 5  # Example DMR ID

input_ids, attention_mask, methylation_tensor, position_tensor, dmr_tensor = prepare_input(
    dna_sequence, methylation_status, position, dmr_id, model.tokenizer
)

# Forward pass
logits = model(input_ids, attention_mask, methylation_tensor, position_tensor, dmr_tensor)
probabilities = torch.softmax(logits, dim=1)
print(f"Probability of tumor: {probabilities[0][1].item():.4f}")
