# LDA experiment

The experiment is executed through the scripts located in the `R/` directory. These scripts should be run sequentially, following the numbering in their filenames (for example, `00_LDA_paths.R` is the first to execute).

Before running these scripts, Wikipedia text dumps are required. We recommend downloading the Wikipedia dumps from [dumps.wikimedia.org](http://dumps.wikimedia.org/enwiki/latest/) and parsing it with [WikiExtractor](https://github.com/attardi/wikiextractor). The resulting files should be placed in the `Data/wikidump` directory. Due to GitHub’s storage limitations, this repository only includes a small sample of documents as an example, downloaded from from [here](https://dumps.wikimedia.org/enwiki/latest/enwiki-latest-pages-articles.xml.bz2), so users need to download the dataset themselves.

The vocabulary used is as in [Hoffman et al. (2010)](https://www.di.ens.fr/~fbach/mdhnips2010.pdf). It is derived from the top 10,000 words in Project Gutenberg texts, excluding all words with fewer than three characters. This results in a vocabulary size of approximately 8,000 words. The vocabulary is provided in `Data/wiki.vocab`.

## Reference

M. Hoffman, F. Bach, and D. Blei. *Online learning for latent Dirichlet allocation.* In J. Lafferty, C. Williams, J. Shawe-Taylor, R. Zemel, and A. Culotta (eds.), Advances in Neural Information Processing Systems, volume 23. Curran Associates, Inc., 2010.

