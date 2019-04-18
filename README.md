# ID3 Decision Tree Classifier with ga4gh_server

An ID3 Decision Tree Classifier that is used to interract with the ga4gh_server (https://github.com/CanDIG/ga4gh-server). The classifier is implemented in such a way so that it allows differential privacy to protect personal health information (PHI).

There is also a ConfusionMatrix class that extends the ID3 class which is used to determine the accuracy of the decision tree.

Note: There is also a local implementation to be used with VCFs installed locally

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

### Prerequisites

What things you need to install the software and how to install them:
```
python3/pip (https://www.python.org/downloads/)
virtualenv (https://virtualenv.pypa.io/en/stable/installation/)
git (https://git-scm.com/)
ga4gh-server (https://github.com/CanDIG/ga4gh-server)
```

### Installing ID3

A step by step series of examples that tell you how to get a development env running

1. Clone the repo locally

```
git clone https://github.com/matthewwong525/ID3_project_sickkids.git
```

2. Setup virtualenv and install dependencies

```
cd ID3_project_sickkids
virtualenv id3_env
source id3_env/bin/activate
pip3 install -r requirements.txt
```

3. Downloading 1000 genomes VCF files to digest (NEEDED FOR LOCAL API)

Download this entire directory (http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/) into a `release` directory in the repo
or run the following commands below

```
wget -m ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ -nd -P release -l 1
```
4. Test if everything works!

Run script below to see if everything works!

```
python main.py
```

End with an example of getting some data out of the system or using it for a little demo

## Updating `config.json`

`config.json` has been already set up so you can run the classifier without any changes to this file. If you want to make some changes to your query, you can simply change some values within this file. Below are the explanations of the attributes within this file.

```
`variant_ranges` : The ranges of variants you want to inspect in the ID3 classifier. A variant range can contain more than one variant.
`ga4gh_server_url` : The url that points to the ga4gh_server
`ga4gh_server_dataset_id` : The id of the dataset you want to query from. It is associated with the ga4gh_server
`user_mapping_path` : Path to the `.ped` file that maps individual ids to ancestries
`chr_paths` : Path to the `.vcf` chromosome files from the 1000 genomes project
```

### Installing and Starting candig_server

Additionally, you can host the ga4gh_server locally if you want to test the capabilities

```
# install OS packages

## for Ubuntu/Debian
sudo apt-get install python-dev python-virtualenv zlib1g-dev libxslt1-dev libffi-dev libssl-dev
export LIBRARY_SEARCH_PATHS=$LIBRARY_SEARCH_PATHS:/usr/local/opt/openssl/lib/

## for Fedora 22+ (current)
sudo dnf install python-devel python-virtualenv zlib-devel libxslt-devel openssl-devel

# setup env
cd ga4gh_server
virtualenv test_server
cd test_server
source bin/activate

# install packages
pip install -U git+https://github.com/CanDIG/candig-schemas.git@develop#egg=ga4gh_schemas
pip install -U git+https://github.com/CanDIG/candig-client.git@authz#egg=ga4gh_client
pip install -U git+https://github.com/CanDIG/candig-server.git@master#egg=candig_server
pip install -U git+https://github.com/CanDIG/candig-ingest@master#egg=candig_ingest
pip install -U git+https://github.com/CanDIG/PROFYLE_ingest.git@develop#egg=PROFYLE_ingest

# setup initial peers
mkdir -p ga4gh/server/templates
touch ga4gh/server/templates/initial_peers.txt

# ingest data and make the repo
mkdir ga4gh-example-data
cd ga4gh-example-data

# init db
ga4gh_repo init registry.db

# create dataset
ga4gh_repo add-dataset registry.db 1kgenome \
    --description "Variants from the 1000 Genomes project and GENCODE genes annotations"

# ingest patient metadata
PROFYLE_ingest registry.db 1kgenome ../../1kgenome_metadata.json

# add reference set
cd ../..
wget ftp://ftp.1000genomes.ebi.ac.uk//vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
cd test_server/ga4gh-example-data

ga4gh_repo add-referenceset registry.db ../../hs37d5.fa.gz \
  -d "NCBI37 assembly of the human genome" --name GRCh37-lite \
  --sourceUri "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz"

# add variant sets (100 variantsets to ingest)
bash ../../1kgenome_ingest.sh

# launch server
cd ..
ga4gh_server --host 127.0.0.1 --port 8000 -c NoAuth

```

## Examples


### ID3 Classifier Example (with 1000 genomes vcf files)

Below is an exanple that uses the 1000 genomes vcf files. If you want to connect to the ga4gh_server, simply change the `local` flag to `False`

This example can be found in main.py

```
from ID3_Class import ID3

# Creates ID3 object with a filepath to the config.json
id3_obj = ID3('config.json', local=True)

# prints the ID3 tree as `tree.png`
id3_obj.print_tree('tree')

# prints the list of all the variant names
print(id3_obj.api.variant_name_list)

# predicts ancestry of the person with the variant `22:50121766:50121767` and no other variant in `variant_name_list`
id3_obj.predict(['22:50121766:50121767'])
```

### Confusion Matrix Example

Confusion Matrix object extends the ID3 object, so you can use functions like `predict`. The Confusion matrix is used to describe performance of the model.

```
from ConfusionMatrix import ConfusionMatrix

# Creates ConfusionMatrix object
c = ConfusionMatrix()

# A few functions from the class that can be used
c.print_matrix()
print c.accuracy()
print c.misclass_rate()
print c.true_ancestry_rate('GBR')
print c.false_ancestry_rate('GBR')
print c.prevalance('GBR')
print c.precision('GBR')
```



