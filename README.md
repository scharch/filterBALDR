# filterBALDR
This is a simple script for removing background noise from
BALDR output before downstream processing.

System Requirements:
* Perl 5 with Pod::Usage and Statistics::Basic modules
   (tested on Perl 5.24.0, Pod::Usage 1.68, 
   		Statistics::Basic 1.6611)

Installation:
< none >

Instructions for use:
* Generate CSV sample sheet with one line perl cell:
	- column 3 is the flow cell ID
	- column 4 is the lane
	- column 5 is the index/sequenceing barcode
	- column 7 is the human readable cell ID
* Run BALDR, with output saved to a path that includes
   <flowcell>/<lane>
* Run BALDR summary scripts
* `perl filterBALDR.pl samplesheet.csv`

For further usage details:
`perl filterBALDR.pl -h`

Expected output:
* Fasta file (`filtered.fa`) of high quality contigs

Run time:
* Typically <1 second
