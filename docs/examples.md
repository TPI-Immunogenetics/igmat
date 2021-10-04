# Examples #

IgMAT can be used to annotate sequences with the default IMGT dataset, or specific datasets can be created:

## Annotating a sequence ##
To annotate an input sequence with the default IMGT database:

	igmat run -i DVQLVESGGGSVQAGGSLRLSCAVSGSTYSPCTTGWYRQAPGKEREWVSSISSPGTIYYQDSVKGRFTISRDNAKNTVYLQMNSLQREDTGMYYCQIQCGVRSIREYWGQGTQVTVSSHHHHHH -v

It will print a representation of the annotated sequence:

```
Query 'Input sequence':

 HMMR match: Vicugna+pacos_H [eValue: 2.000e-55]

   0 -----------FR1------------|----CDR1----|-------FR2-------|---CDR2---|-  70
   0 DVQLVESGG-GSVQAGGSLRLSCAVS|GSTY----SPCT|TGWYRQAPGKEREWVSS|ISSP---GTI|Y  70

  70 -----------------FR3------------------|-----CDR3----|----FR4----| 135
  70 YQDSVK-GRFTISRDNAKNTVYLQMNSLQREDTGMYYC|QI-QCGVRSIREY|WGQGTQVTVSS| 135
```

## Creating a custom dataset ##
To create a custom dataset, input files with sequences from V and J regions are needed, where the name of the files must adhere the following format: 
	
		<name>_<chain><type>
Where `name` is a species or sample name, and `type` is either *V* or *J*.
An example of input sequences for generating custom datasets is distributed in the `test/build` folder:

	igmat build -i ./test/build -n cattle
The command will generate a `cattle` dataset that will be available to be used with IgMAT with the command:
	
	igmat run -i DVQLVESGGGSVQAGGSLRLSCAVSGSTYSPCTTGWYRQAPGKEREWVSSISSPGTIYYQDSVKGRFTISRDNAKNTVYLQMNSLQREDTGMYYCQIQCGVRSIREYWGQGTQVTVSSHHHHHH -m cattle