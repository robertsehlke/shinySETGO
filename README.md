# shinySETGO
Shiny tool, Simple Enrichment Testing for Gene Ontology.

This is a side-project mainly intended to tinker with the shiny framework,
and refresh my knowledge of mySQL databases (and its integration with R).

The result is a shiny web app that performs GO enrichment analysis of gene sets of
interest, currently in two organisms: Drosophila and C. elegans. The user can
provide genes as entrez identifiers, flybase/wormbase ID, or gene names. It is also
possible to specify a custom background set of genes, e.g. in cases where the method of
collection already biases the functional content (tissue- or cell-specific data, etc.).

To cut down on redundancy, the app (through my SETHRO script, see repo here) calculates
clusters of GO terms that differ in their annotated genes of interest by at most 5 genes
(this will be adjustable later). The most specific term (usually the smallest) is then
selected as a representative term for each cluster. Testing for significant enrichment
is performed using a Fisher test, after which p-values of primary terms are adjusted for
multiplicity using the Benjamini-Hochberg method. 

This a simpler approach than existing tools like DAVID and topGO that implement clustering
or elimination of redundant terms. The benefit in my own projects was to keep the
workflow completely in R (unlike DAVID) and have a simple and transparent reduction of
redundancy (unlike topGOs elim approach). Ultimately it's a small tweak to a known
procedure that proved to be helpful to myself, as well as a good learning experience.


## Note about annotations:
Annotations for Drosophila and C.elegans were taken from the following bioconductor 
annotation packages:
https://bioconductor.org/packages/release/data/annotation/html/org.Dm.eg.db.html
* Carlson M (2017). org.Dm.eg.db: Genome wide annotation for Fly. R package version 3.5.0.
https://bioconductor.org/packages/release/data/annotation/html/org.Ce.eg.db.html
* Carlson M (2017). org.Ce.eg.db: Genome wide annotation for Worm. R package version 3.5.0.

Additional mappings were taken from BioMart (http://www.biomart.org).

The data necessary for this tool were then transferred into a custom database scheme.
