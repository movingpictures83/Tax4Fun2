# Tax4Fun2
# Language: R
# Input: TXT (keyword, value pairs)
# Output: DIR (for results)
# Tested with: PluMA 1.1, R 4.0.0
# Dependency: Tax4Fun2 

Note: Tax4Fun2 is available unofficially, open source at: https://github.com/ZihaoShu/Tax4Fun2

This plugin was partially built using a modified version of that package.

The plugin takes a TXT file of tab-delimited keyword-value pairs:

fasta: FASTA file of OTU representative sequences.
otutable: Tab-delimited OTU table (OTUs are rows, samples are columns)
pathtoreference: Path to Tax4Fun2 reference database (you'll need to install this using the above URL)

The plugin outputs both functional and pathway prediction TXT files in the user-specified output directory.
