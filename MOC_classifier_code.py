#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# =====================================================
# Proteomic Pan-Cancer Classifier
# Multi-Model Framework (Aligned to Fig 2)
# =====================================================

# This notebook implements multiple modelling frameworks:
# 1. Proof-of-principle (internal validation)
# 2. Expanded training + cohort-based validation
# 3. External validation (Qian)
# 4. Prediction + interpretation layers


# In[ ]:


# Input:
#   combat_Omicodf5625.xlsx as derived from initial workflow and detailed in accompanying R GitHub repository 
#
# Feature panel:
#   Supplementary Table 4 (37 proteins)
#
# Outputs:
#   Figure 2B-E
#   Supplementary Figure 3A-D
#
# Models:
#   Model 1: Internal validation
#   Model 2: Expanded training
#   Model 3: Omico holdout validation
#   External validation: Qian cohort

# The classifier uses a predefined 37-protein diagnostic panel
# described in Supplementary Table 4.
#
# This panel was finalised prior to classifier training and is
# not generated within this script.


# In[ ]:


import pandas as pd
import numpy as np 
from sklearn import preprocessing
from sklearn.model_selection import train_test_split
from sklearn.feature_selection import SelectKBest
from sklearn.feature_selection import f_regression
import seaborn as sns
from sklearn import datasets, metrics, model_selection, svm
from sklearn.preprocessing import OneHotEncoder 
from sklearn.metrics import top_k_accuracy_score
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import cross_val_score
import matplotlib.pyplot as plt
get_ipython().run_line_magic('matplotlib', 'inline')
from sklearn.metrics import confusion_matrix
from sklearn.preprocessing import label_binarize
from sklearn.preprocessing import LabelEncoder
from imblearn.over_sampling import RandomOverSampler
from sklearn.model_selection import GroupKFold, StratifiedKFold
from sklearn.feature_selection import RFECV
from sklearn.ensemble import RandomForestClassifier, VotingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.pipeline import Pipeline
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import accuracy_score, roc_curve, auc
from sklearn.feature_selection import RFECV
from sklearn.svm import SVC
from xgboost import XGBClassifier
from sklearn.pipeline import Pipeline
from sklearn.metrics import balanced_accuracy_score, roc_auc_score, classification_report, f1_score, make_scorer
from sklearn.metrics import precision_recall_curve, average_precision_score
from itertools import cycle
from sklearn.multiclass import OneVsRestClassifier


# In[ ]:


# Read the excel file with additional parameters
df1 = pd.read_excel('combat_Omicodf5625.xlsx', header=0)


# In[ ]:


# Filter for relevant diagnoses
df2 = df1[~df1['Label'].isin([
    'Appendix_adenocarcinoma',
    'Small_bowel_adenocarcinoma',
    'Bladder_mucinous'
])]


# In[ ]:


# Fixed 37-protein diagnostic panel
# Derived during feature selection analyses described in
# Supplementary Table 4.
# No feature selection is performed in this notebook.

selected_features = ['Q12864','Q9UJ72','P20142','Q9UPW6','Q99795','P56856','P01275','Q8NBJ4','Q86XK7','Q16790',
                     'P19075','P06731','P25815','P19544','O94819','Q9HCB6','Q8IUD2','P04626','Q2TBA0','Q9BYZ2',
                     'Q9BSA9','Q9BYP7','Q6UY09','A1Z1Q3','Q10472','Q9UPW5','Q16625','P28330','O60513','P36222',
                     'P04844','A8K7I4','P54687', 
                    'Q9UBC9', 'P07098', 'P08729', 'P35900']


# In[ ]:


# Remove rows where sample_cohort starts with "EC" (Omico sample set)
df4 = df2[~df2['sample_cohort'].str.startswith('EC', na=False)]


# In[ ]:


# Assign variables
X = df4.iloc[:, 0:10416]  # Features in columns 0:10415
y = df4['Label']  # Target (diagnosis)
groups = df4['sample_cohort']  # Cohort information

# Encode labels
le = LabelEncoder()
y_encoded = le.fit_transform(y)

# Check class counts
class_counts = y.value_counts()
class_counts


# In[ ]:


X_selected = X[selected_features]

# Create a stratification variable
stratify_var = pd.Series(list(zip(y)))

# Split the data with stratification
X_train, X_test, y_train, y_test = train_test_split(X_selected, y_encoded, test_size=0.2, random_state=42, stratify=stratify_var)



# In[ ]:


# Fit Model 1 - internal model 
svc_clf = SVC(probability=True, random_state=42)
svc_param_grid = {
    'C': [1, 10],
    'gamma': ['scale', 0.01]
}

svc_grid = GridSearchCV(svc_clf, svc_param_grid, cv=3, n_jobs=-1, verbose=2)
svc_grid.fit(X_train, y_train)

# Define the original labels
labels = ['Colorectal adenocarcinoma', 'High grade serous ovary', 
          'Mucinous ovarian', 'Pancreas adenocarcinoma', 'Gastro-oesophageal adenocarcinoma']


print("Best SVC Params:", svc_grid.best_params_)
y_pred_svc = svc_grid.predict(X_test)
print("SVC Accuracy:", accuracy_score(y_test, y_pred_svc))
print("SVC Classification Report:\n", classification_report(y_test, y_pred_svc, target_names=labels))


# In[ ]:


# Supplementary Fig 3A - confusion  matrix for internal model 
from sklearn.metrics import ConfusionMatrixDisplay, confusion_matrix

labels = ['CRC', 'HGSOC', 'MOC', 'PDAC', 'GEJ']
cm = confusion_matrix(y_test, y_pred_svc)

# Convert to row-wise proportions
row_totals = cm.sum(axis=1, keepdims=True)
cm_decimal = cm / row_totals
cm_formatted = np.where(cm == 0, "0", np.round(cm_decimal, 2).astype(str))  # handle zeros

# Create display object
# Normalize to row-wise proportions
row_totals = cm.sum(axis=1, keepdims=True)
cm_normalized = cm / row_totals

# Create display object with normalized values
disp = ConfusionMatrixDisplay(confusion_matrix=cm_normalized, display_labels=labels)
# Plot
fig, ax = plt.subplots(figsize=(6, 6))
disp.plot(ax=ax, xticks_rotation='vertical', cmap='Blues')

# Replace text with formatted decimals
for i in range(cm.shape[0]):
    for j in range(cm.shape[1]):
        disp.text_[i, j].set_text(cm_formatted[i, j])
        disp.text_[i, j].set_fontsize(14)
        disp.text_[i, j].set_fontweight('bold')

# Axis labels
ax.set_xlabel('Predicted label', fontsize=14, weight='bold')
ax.set_ylabel('True label', fontsize=14, weight='bold')
ax.tick_params(axis='both', labelsize=14)

plt.tight_layout()
plt.savefig("confusion_matrix_small.jpg", dpi=300, format='jpg', bbox_inches='tight')
plt.show()


# In[ ]:


# Supplementary Fig 3b - ROC-AUC curves for internal model 
# Assume you have:
# X_train_selected, y_train, X_test_selected, y_test, labels

# 1. Binarize the output for ROC
y_test_bin = label_binarize(y_test, classes=np.arange(len(labels)))
n_classes = y_test_bin.shape[1]

# 3. Predict probabilities
y_score = svc_grid.predict_proba(X_test)

# 4. Plot ROC curve for each class
plt.figure(figsize=(8, 6))
for i in range(n_classes):
    fpr, tpr, _ = roc_curve(y_test_bin[:, i], y_score[:, i])
    roc_auc = auc(fpr, tpr)
    plt.plot(fpr, tpr, lw=2, label=f'Class {labels[i]} (AUC = {roc_auc:.2f})')
plt.plot([0, 1], [0, 1], 'k--', lw=2)
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('One-vs-Rest ROC Curves for SVC')
plt.legend(loc="lower right")
plt.show()


# In[ ]:


# Supplementary Fig 3c - Precision recall curves for internal model 

# Binarize labels
y_test_bin = label_binarize(y_test, classes=np.arange(len(labels)))
n_classes = y_test_bin.shape[1]

# Predict probabilities
y_score = best_svc.predict_proba(X_test)

# Compute PR curves and average precision
precision, recall, average_precision = dict(), dict(), dict()
for i in range(n_classes):
    precision[i], recall[i], _ = precision_recall_curve(y_test_bin[:, i], y_score[:, i])
    average_precision[i] = average_precision_score(y_test_bin[:, i], y_score[:, i])

# Micro-average
precision["micro"], recall["micro"], _ = precision_recall_curve(y_test_bin.ravel(), y_score.ravel())
average_precision["micro"] = average_precision_score(y_test_bin, y_score, average="micro")

# Plot
fig, ax = plt.subplots(figsize=(8, 6))
colors = cycle(["aqua", "darkorange", "cornflowerblue", "green", "red"])

# Micro-average curve
ax.plot(
    recall["micro"],
    precision["micro"],
    label=f"Micro-average (AP = {average_precision['micro']:.2f})",
    color="deeppink",
    linestyle=":",
    linewidth=4,
)

# Per-class curves
for i, color in zip(range(n_classes), colors):
    ax.plot(
        recall[i],
        precision[i],
        color=color,
        lw=2,
        label=f"{labels[i]} (AP = {average_precision[i]:.2f})"
    )

# Axis styling
ax.set_xlabel("Recall", fontsize=14, weight='bold')
ax.set_ylabel("Precision", fontsize=14, weight='bold')
ax.tick_params(axis='both', labelsize=14)
ax.set_facecolor('white')

# Black axis lines only
ax.spines['bottom'].set_color('black')
ax.spines['left'].set_color('black')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Legend outside lower right
legend = ax.legend(
    loc="lower left",
    fontsize=14,
    frameon=False
)

# Export
plt.tight_layout()
plt.savefig("precision_recall_curve_small.jpg", dpi=300, format='jpg', bbox_inches='tight', bbox_extra_artists=[legend])
plt.show()


# In[ ]:


# CV score for internal model
from sklearn.model_selection import cross_val_score

best_svc = svc_grid.best_estimator_
cv_scores = cross_val_score(best_svc, X_train, y_train, cv=5, scoring='accuracy')

print("5-fold CV accuracies:", cv_scores)
print("Mean 5-fold CV accuracy: {:.3f}".format(cv_scores.mean()))
print("Standard deviation: {:.3f}".format(cv_scores.std()))


# In[ ]:


# =====================================================
# MODEL 2 — EXPANDED TRAINING MODEL (FIG 2b–d)
# Cohort 1 + 2 + 50% Omico for training
# Remaining 50% Omico used as independent test set
# =====================================================
# MODEL 2 — FINAL PERFORMANCE MODEL (FIG 2B-D)
# --- Define dataset ---
X = df2.iloc[:, 0:10416]
y = df2['Label']

le = LabelEncoder()
y_encoded = le.fit_transform(y)

# --- Identify Omico cohort ---
mask_ec = df2['sample_cohort'].str.startswith('EC', na=False)

# Split only NON-MOC EC samples (as per your logic)
mask_mucinous = (df2['Label'] == 'Mucinous_ovarian') & mask_ec
mask_ec_nonmucinous = mask_ec & ~mask_mucinous

# --- Split EC (Omico) into train/test ---
from sklearn.model_selection import train_test_split

idx_ec = np.where(mask_ec_nonmucinous)[0]

train_idx, test_idx = train_test_split(
    idx_ec,
    test_size=0.5,
    random_state=42,
    stratify=y_encoded[idx_ec]
)

# --- Construct final masks ---
mask_non_ec = ~mask_ec

final_train_mask = mask_non_ec.to_numpy().copy()
final_train_mask[train_idx] = True

final_test_mask = np.zeros_like(mask_ec, dtype=bool)
final_test_mask[test_idx] = True
final_test_mask[mask_mucinous.to_numpy()] = True  # include MOC in test

# --- Apply feature selection ---
X_selected = X[selected_features]

X_train = X_selected[final_train_mask]
y_train = y_encoded[final_train_mask]

X_test = X_selected[final_test_mask]
y_test = y_encoded[final_test_mask]

# --- Train model ---
from sklearn.svm import SVC

svc = SVC(probability=True, random_state=42)

param_grid = {
    'C': [1, 10],
    'gamma': ['scale', 0.01],
    'kernel': ['rbf']
}

grid = GridSearchCV(
    svc,
    param_grid,
    cv=5,
    n_jobs=-1,
    scoring='accuracy'
)

grid.fit(X_train, y_train)

best_svc = grid.best_estimator_

# --- Evaluate ---

y_pred = best_svc.predict(X_test)

print("Test accuracy:", accuracy_score(y_test, y_pred))
print(classification_report(y_test, y_pred))

# =====================================================
# FIG 2b — CONFUSION MATRIX
# =====================================================

from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay

labels = ['CRC', 'HGSOC', 'MOC', 'PDAC', 'GEJ']

cm = confusion_matrix(y_test, y_pred)
cm_norm = cm / cm.sum(axis=1, keepdims=True)

disp = ConfusionMatrixDisplay(confusion_matrix=cm_norm, display_labels=labels)

fig, ax = plt.subplots(figsize=(6,6))
disp.plot(ax=ax, cmap='Blues', colorbar=False)
plt.tight_layout()
plt.show()

# =====================================================
# FIG 2c — ROC CURVES
# =====================================================
y_test_bin = label_binarize(y_test, classes=np.arange(len(labels)))
y_score = best_svc.predict_proba(X_test)

fpr, tpr, roc_auc = {}, {}, {}

for i in range(len(labels)):
    fpr[i], tpr[i], _ = roc_curve(y_test_bin[:, i], y_score[:, i])
    roc_auc[i] = auc(fpr[i], tpr[i])

# micro-average
fpr["micro"], tpr["micro"], _ = roc_curve(y_test_bin.ravel(), y_score.ravel())
roc_auc["micro"] = auc(fpr["micro"], tpr["micro"])

# plot
plt.figure(figsize=(8,6))
for i in range(len(labels)):
    plt.plot(fpr[i], tpr[i], label=f"{labels[i]} (AUC={roc_auc[i]:.2f})")
plt.plot(fpr["micro"], tpr["micro"], linestyle="--", label=f"Micro AUC={roc_auc['micro']:.2f}")
plt.xlabel("FPR"); plt.ylabel("TPR")
plt.legend()
plt.show()

# =====================================================
# FIG 2d — PRECISION-RECALL
# =====================================================

from sklearn.metrics import precision_recall_curve, average_precision_score

precision, recall, ap = {}, {}, {}

for i in range(len(labels)):
    precision[i], recall[i], _ = precision_recall_curve(y_test_bin[:, i], y_score[:, i])
    ap[i] = average_precision_score(y_test_bin[:, i], y_score[:, i])

plt.figure(figsize=(8,6))
for i in range(len(labels)):
    plt.plot(recall[i], precision[i], label=f"{labels[i]} (AP={ap[i]:.2f})")
plt.xlabel("Recall"); plt.ylabel("Precision")
plt.legend()
plt.show()


# In[ ]:


# SHAP Analysis - derived from Model 2 - expanded training dataset
# This analysis used to create Figure 2E
import shap

# 1. Define a prediction function for your trained model
def predict_fn(X):
    return best_svc.predict_proba(X)

# 2. Use a sample of your training data as the background (for performance, use a small random subset)
background = X_train.iloc[np.random.choice(X_train.shape[0], 100, replace=False)]

# 3. Create the SHAP KernelExplainer
explainer = shap.KernelExplainer(predict_fn, background)

# 4. Compute SHAP values for your test set (or a subset, for speed)
shap_values = explainer.shap_values(X_test[:174])  # Use a manageable subset for speed


# In[ ]:


# Transpose SHAP values to (n_classes, n_samples, n_features)
shap_values_per_class = np.transpose(shap_values, (2, 0, 1))  # shape: (n_classes, n_samples, n_features)

# Loop through each class
for i, class_shap in enumerate(shap_values_per_class):
    # Convert to DataFrame
    shap_df = pd.DataFrame(class_shap, columns=X_test.columns[:class_shap.shape[1]])
    
    # Add sample IDs and predicted class
    shap_df['Sample_ID'] = X_test.index[:class_shap.shape[0]].values
    shap_df['Predicted_Class'] = best_svc.predict(X_test[:class_shap.shape[0]])
    
    # Reorder columns
    cols = ['Sample_ID', 'Predicted_Class'] + list(X_test.columns[:class_shap.shape[1]])
    shap_df = shap_df[cols]
    
    # Export to Excel
    shap_df.to_excel(f"shap_values_class_{i}.xlsx", index=False)


# In[ ]:


# Transpose to (n_classes, n_samples, n_features)
shap_values_per_class = np.transpose(shap_values, (2, 0, 1))

# Create a dictionary to hold mean SHAP values per class
shap_summary = {}

for i, class_shap in enumerate(shap_values_per_class):
    # Compute mean(|SHAP|) across samples for each feature
    shap_summary[f'Class_{i}'] = np.abs(class_shap).mean(axis=0)

# Convert to DataFrame: rows = features, columns = classes
shap_df = pd.DataFrame(shap_summary, index=X_test.columns)

# Reset index to make 'Feature' a column
shap_df.reset_index(inplace=True)
shap_df.rename(columns={'index': 'Feature'}, inplace=True)

# Export to Excel
shap_df.to_excel("shap_feature_summary.xlsx", index=False)


# In[ ]:


# Convert to (n_classes, n_samples, n_features)
shap_values_per_class = np.transpose(shap_values, (2, 0, 1))
# Now shap_values_per_class[0] is (100, 37) for class 0, etc.


# In[ ]:


#from UniProtMapper import ProtMapper

mapper = ProtMapper()

# Replace with your list of UniProt IDs
uniprot_ids = ['Q12864','Q9UJ72','P20142','Q9UPW6','Q99795','P56856','P01275','Q8NBJ4','Q86XK7','Q16790',
                     'P19075','P06731','P25815','P19544','O94819','Q9HCB6','Q8IUD2','P04626','Q2TBA0','Q9BYZ2',
                     'Q9BSA9','Q9BYP7','Q6UY09','A1Z1Q3','Q10472','Q9UPW5','Q16625','P28330','O60513','P36222',
                     'P04844','A8K7I4','P54687', 
                    'Q9UBC9', 'P07098', 'P08729', 'P35900']

# Map from UniProt to Gene Names (primary symbols)
result, failed = mapper.get(
    ids=uniprot_ids,
    from_db="UniProtKB_AC-ID",
    to_db="Gene_Name"
)

print(result)  # pandas DataFrame with mapping results
print(failed)  # list of IDs that couldn't be mapped


# In[ ]:


# SHAP summary plots
# For class 0
shap.summary_plot(shap_values_per_class[0], X_test[:100], feature_names=gene_symbols)

# For class 1
shap.summary_plot(shap_values_per_class[1], X_test[:100], feature_names=gene_symbols)

# For class 2
shap.summary_plot(shap_values_per_class[2], X_test[:100], feature_names=gene_symbols)

# For class 3
shap.summary_plot(shap_values_per_class[3], X_test[:100], feature_names=gene_symbols)

# For class 4
shap.summary_plot(shap_values_per_class[4], X_test[:100], feature_names=gene_symbols)


# In[ ]:


# =====================================================
# MODEL 3 — OMICO FULL HOLDOUT (SUPPLEMENTARY FIG 3d)
# Train on cohorts 1 + 2
# Test on ALL Omico samples (cohort 3)
# =====================================================

# IMPORTANT:
# Cohort 3 (Omico) was NOT used for feature selection.
# This dataset is used as an independent validation cohort
# to assess generalisability and mitigate potential overfitting.

# MODEL 3 — INDEPENDENT VALIDATION (SUPP FIG 3d)

# --- Define dataset ---
X = df2.iloc[:, 0:10416]
y = df2['Label']

le = LabelEncoder()
y_encoded = le.fit_transform(y)

# --- Split cohorts ---
mask_ec = df2['sample_cohort'].str.startswith('EC', na=False)

X_selected = X[selected_features]

X_train = X_selected[~mask_ec]
y_train = y_encoded[~mask_ec]

X_test = X_selected[mask_ec]
y_test = y_encoded[mask_ec]

# --- Train ---
svc = SVC(probability=True, random_state=42)

grid = GridSearchCV(
    svc,
    {'C':[1,10], 'gamma':['scale',0.01]},
    cv=5,
    n_jobs=-1
)

grid.fit(X_train, y_train)
best_svc = grid.best_estimator_

# --- Evaluate ---
y_pred = best_svc.predict(X_test)

print("Omico holdout accuracy:", accuracy_score(y_test, y_pred))
print(classification_report(y_test, y_pred))

# =====================================================
# SUP FIG 3d — CONFUSION MATRIX
# =====================================================

labels = ['CRC', 'HGSOC', 'MOC', 'PDAC', 'GEJ']

cm = confusion_matrix(y_test, y_pred)
cm_norm = cm / cm.sum(axis=1, keepdims=True)

disp = ConfusionMatrixDisplay(confusion_matrix=cm_norm, display_labels=labels)

fig, ax = plt.subplots(figsize=(6,6))
disp.plot(ax=ax, cmap='Blues')
plt.title("Omico full holdout")
plt.show()
``


# In[ ]:


# =====================================================
# Merge PDAC + GEJ → Upper GI
# Done due to limited sample numbers (<20 in training data)
# =====================================================

df2['Label_merged'] = df2['Label'].replace({
    'Pancreas_adenocarcinoma': 'Upper GI',
    'Gastro-oesophageal_adenocarcinoma': 'Upper GI'
})

# Note:
# HGSOC not included in this model due to insufficient representation
# in the training dataset to support robust classification
# Confusion matrix plot as per above code


# In[ ]:


# =====================================================
# CLINICAL APPLICATION MODEL — METASTATIC OVARIAN SAMPLES
# Model retrained excluding selected test samples
# NOT USED FOR PERFORMANCE METRICS
# =====================================================

# Define metastatic / clinical test samples
test_samples = ['03NKT', '03NNL', '040YB', '03QYU',
                '03NNH', '03QXY', '03XM7', '03XLV',
                '04JO2', '040Y6']

# Split dataset
test_set = df3[df3['Sample_ID'].isin(test_samples)]
train_set = df3[~df3['Sample_ID'].isin(test_samples)]

# --- Feature selection ---
X_train = train_set[selected_features]
y_train = train_set['Label']

X_test = test_set[selected_features]
y_test = test_set['Label']

# --- Encode labels ---
le = LabelEncoder()
y_train_enc = le.fit_transform(y_train)
y_test_enc = le.transform(y_test)

# --- Train model (retrained classifier) ---
from sklearn.svm import SVC
from sklearn.model_selection import GridSearchCV

svc = SVC(probability=True, random_state=42)

param_grid = {
    'C': [1, 10],
    'gamma': ['scale', 0.01],
    'kernel': ['rbf']
}

grid = GridSearchCV(
    svc,
    param_grid,
    cv=5,
    n_jobs=-1,
    scoring='accuracy'
)

grid.fit(X_train, y_train_enc)

best_svc_clinical = grid.best_estimator_

# --- Predict on held-out clinical samples ---
y_pred = best_svc_clinical.predict(X_test)
y_proba = best_svc_clinical.predict_proba(X_test)

# --- Output ---
results = pd.DataFrame({
    'Sample_ID': test_set['Sample_ID'],
    'True_Label': y_test,
    'Predicted_Label': le.inverse_transform(y_pred)
})

print(results)


# In[ ]:


# Apply final model to appendiceal case 
sample = df1[df1['Sample_ID'] == '040Y6']
X_sample = sample[selected_features]

predicted_class = best_svc_clinical.predict(X_sample)
predicted_proba = best_svc_clinical.predict_proba(X_sample)


# In[ ]:


# =====================================================
# EXTERNAL VALIDATION — QIAN COHORT (RETRAINED MODEL)
# Qian excluded from training, used as independent test set
# =====================================================
# External cohort independent of feature selection and internal training
# Used to assess generalisability of classifier

# --- Load Qian dataset ---
df_qian = pd.read_csv("Qian_combat_MOC_dataset.csv")

# --- Combine datasets (excluding Qian from training) ---
train_set = df2  # your full dataset (cohort 1+2+3)
test_set  = df_qian

# --- Feature selection ---
X_train = train_set[selected_features]
y_train = train_set['Label']

X_test = test_set[selected_features]
y_test = test_set['Label']

# --- Encode labels ---
le = LabelEncoder()
y_train_enc = le.fit_transform(y_train)
y_test_enc = le.transform(y_test)

# --- Train model (retrained classifier) ---
from sklearn.svm import SVC

svc = SVC(probability=True, random_state=42)

param_grid = {
    'C': [1, 10],
    'gamma': ['scale', 0.01],
    'kernel': ['rbf']
}

grid = GridSearchCV(
    svc,
    param_grid,
    cv=5,
    n_jobs=-1,
    scoring='accuracy'
)

grid.fit(X_train, y_train_enc)

best_svc_qian = grid.best_estimator_

# --- Predict ---
y_pred = best_svc_qian.predict(X_test)
y_proba = best_svc_qian.predict_proba(X_test)

# --- Evaluate ---
print("Qian Classification Report:")
print(classification_report(y_test_enc, y_pred))

# --- MOC-specific metrics ---
labels = le.classes_
moc_index = list(labels).index('Mucinous_ovarian')

y_true_moc = (y_test_enc == moc_index).astype(int)
y_pred_moc = (y_pred == moc_index).astype(int)

print("MOC Recall:", recall_score(y_true_moc, y_pred_moc))
print("MOC Precision:", precision_score(y_true_moc, y_pred_moc))

