# BOE: BiOnt Enhanced (DrugProt corpus)

The motivation behind this project was to get better results, on the [DrugProt corpus](https://zenodo.org/record/5042151) (drug and chemical-protein interactions), through the adaptation of the [BiOnt](https://github.com/lasigeBioTM/BiOnt) system.

## Installation

Pursue the steps of requirements.txt in the folders: **bin**; **data**; **drugprot-evaluation-library**.

## Enhancements

- Ontology mapping and Ancestor linking are multiprocessed;
- Added faster sentence segmenter (pysbd);
- GO and ChEBI entities ancestors are populated through the .db file and via API;
- Creation of additional senteces with the overlapped entities;
- Improved text tokenization;
- Optimization of some functions;
- Fixed minor bugs.

## Results (Micro-averaged)

Principal variations between the models shown below are as follows:

[**baseline**](https://github.com/lasigeBioTM/biocreativeVII): Model generated using a modified version of the BiOnt system used in the BioCreative VII Track 1 challenge.

**main ents**: The largest term from the previously identified overlapping entities is used, while the rest are discarded.

**extra ents**: Creates new sentences that incorporate one GENE term and one CHEMICAL term mentioned in the original sentence.

| Model | Precision | Recall | F1-Score |
| ------ | :------: | :------: | :------: |
| baseline | **0.618** | 0.097 | 0.168 |
| main ents | 0.24 | 0.293 | 0.264 |
| extra ents | 0.149 | 0.262 | 0.19 |
| *main ents (updated)* | 0.444 | **0.447** | **0.445** |
| *extra ents (updated)* | 0.449 | 0.431 | 0.44 |

Micro-averaged results of the three models predicting the relations of the DrugProt development dataset.

## Preprocessing

- $2: type_of_action (preprocess) 
- $3: preprocessing_method (extra or main) 
- $4: temporary_directory

### Example:

```sh
python src/ontologies_embeddings.py preprocess main temp/

python src/ontologies_embeddings.py preprocess main temp_dev/
```

## Training

- $2: type_of_action (train) 
- $3: model_name 
- $4: channels (words; wordnet; concatenation_ancestors)
- $5: temporary_directory

### Example:

```sh
python src/ontologies_embeddings.py train model_name words wordnet concatenation_ancestors temp/
```

## Predicting

- $2: type_of_action (test)
- $3: model_name 
- $4: channels (words; wordnet; concatenation_ancestors)
- $5: temporary_directory

### Example:

```sh
python src/ontologies_embeddings.py test model_name words wordnet concatenation_ancestors temp_dev/
```

## Evaluating

Coverting the result file into a friendly format to be evaluated

- $2: path to Prediction TXT file
- $3: path to new Prediction TSV file

### Example:

```sh
python scripts/validation.py results/model_name_results.txt results/model_name_results.tsv
```

Using the [DrugProt Evaluation library]( https://github.com/tonifuc3m/drugprot-evaluation-library)

- $2: -g/--gs_path: path to Gold Standard relations TSV file 
- $3: -p/--pred_path: path to Prediction TSV file
- $4: -e/--ent_path: path to Gold Standard entities TSV file 
- $5: --pmids: path to list of relevant PMIDs

### Example:

```sh
python drugprot-evaluation-library/src/main.py -g gs-data/drugprot_development_relations.tsv -p results/model_name_results.tsv -e gs-data/drugprot_development_entities.tsv --pmids gs-data/pmids.txt
```
