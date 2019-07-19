## BADGE: Batch Analysis of Differential Gene Expression

BADGE is an open source tool for determining genes that are differentially expressed.

### Requirements
``` shell
python 3.6.x
numpy
pandas
scipy
matplotlib
```

### Setup

Create new directory for project:
``` shell
mkdir \path\to\dir
cd \path\to\dir
```

Create virtual environment:
``` shell
virtualenv badgevenv --python=python3.6
source ./badgevenv/bin/activate
```

Clone the directory and make it the current directory:
``` shell
git clone https://github.com/avkondepudi/badge.git
cd badge
```

Download required dependencies:
``` shell
pip install -r requirements.txt
```

### Usage

#### Determining Differentially Expressed Genes:
``` shell
usage: genes.py [-h] [-g GRAPH] input [ctypes [ctypes ...]]

Get differentially expressed genes.

positional arguments:
input                 input file
ctypes                cell types

optional arguments:
-h, --help            show this help message and exit
-g GRAPH, --graph GRAPH
whether to graph (default: True)
```
