import torch
from transformers import BertForSequenceClassification
model = BertForSequenceClassification.from_pretrained('Peltarion/dnabert-minilm')

