# BiOnt

## Preprocessing
### $2: type_of_action $3: pair_type $4: preprocess_what $5: input_path $6: temporary_directory

#### python3 src/ontologies_embeddings.py preprocess DRUG-GENE train corpora/drug_gene/train temp/

#### python3 src/ontologies_embeddings.py preprocess DRUG-GENE test corpora/drug_gene/test temp_dev/


## Training
### $2: type_of_action $3: pair_type $4: model_name $6: channels $7: temporary_directory

#### python3 src/ontologies_embeddings.py train DRUG-GENE model words wordnet concatenation_ancestors temp/


## Predicting
### $2: type_of_action $3: pair_type $4: model_name $5: gold_standard $6: channels $7: temporary_directory

#### python3 src/ontologies_embeddings.py test DRUG-GENE model corpora/drug_gene/test/ words wordnet concatenation_ancestors temp_dev/


## Evaluating
### -g/--gs_path: path to Gold Standard relations TSV file -p/--pred_path: path to Prediction TSV file -e/--ent_path: path to Gold Standard entities TSV file --pmids: path to list of relevant PMIDs

#### python3 container/drugprot-evaluation-library-main/src/main.py -g container/biont/gs-data/drugprot_development_relations.tsv -p container/biont/results/model.tsv -e container/biont/gs-data/drugprot_development_entities.tsv --pmids container/biont/gs-data/pmids.txt
