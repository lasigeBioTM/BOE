# BOE: BiOnt Enhanced (DrugProt corpus)

Download the [DrugProt corpus](https://zenodo.org/record/5042151)

The aim of this project was to get better results on the DrugProt (drug and chemical-protein interactions) corpus, through the adaptation of the [BiOnt](https://github.com/lasigeBioTM/BiOnt) system, by enhancing it with improved ontologies from Gene Ontology (GO) and Chemical Entities of Biological Interest (ChEBI), and by reducing significantly the amount of time spent on the preprocessing steps. 

It can now capture more entities, allowing overlapping entities to be recognized by spacy, but unfortunately this step failed to improve the results.

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
| ------ | ------ | ------ | ------ |
| baseline | **0.618** | 0.097 | 0.168 |
| main ents | 0.24 | **0.293** | **0.264** |
| extra ents | 0.149 | 0.262 | 0.19 |

Micro-averaged results of the three models predicting the relations of the DrugProt development dataset.

## Preprocessing

- $2: type_of_action 
- $3: pair_type 
- $4: preprocess_what
- $5: input_path 
- $6: temporary_directory

#### Example:

```sh
python3 src/ontologies_embeddings.py preprocess DRUG-GENE train corpora/drug_gene/train temp/
python3 src/ontologies_embeddings.py preprocess DRUG-GENE test corpora/drug_gene/test temp_dev/
```

## Training

- $2: type_of_action 
- $3: pair_type 
- $4: model_name 
- $6: channels 
- $7: temporary_directory

#### Example:

```sh
python3 src/ontologies_embeddings.py train DRUG-GENE model words wordnet concatenation_ancestors temp/
```

## Predicting

- $2: type_of_action 
- $3: pair_type 
- $4: model_name 
- $5: gold_standard 
- $6: channels 
- $7: temporary_directory

#### Example:

```sh
python3 src/ontologies_embeddings.py test DRUG-GENE model corpora/drug_gene/test/ words wordnet concatenation_ancestors temp_dev/
```

## Evaluating

Coverting the result file into a friendly format to be evaluated

- $2: path to Prediction TXT file
- $3: path to new Prediction TSV file

#### Example:

```sh
python3 scripts/validation.py results/model_drug_gene_results.txt results/model_class_weights_drug_gene_results.tsv
```

Using the [DrugProt Evaluation library]( https://github.com/tonifuc3m/drugprot-evaluation-library)

- $2: -g/--gs_path: path to Gold Standard relations TSV file 
- $3: -p/--pred_path: path to Prediction TSV file
- $4: -e/--ent_path: path to Gold Standard entities TSV file 
- $5: --pmids: path to list of relevant PMIDs

#### Example:

```sh
python3 drugprot-evaluation-library-main/src/main.py -g gs-data/drugprot_development_relations.tsv -p results/model.tsv -e gs-data/drugprot_development_entities.tsv --pmids container/biont/gs-data/pmids.txt
```
