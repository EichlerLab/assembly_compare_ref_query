# assembly_compare_ref_query
Assess the QV of an assembly

This pipeline finds differences between assemblies using bwa, then uses freebayes calls to identify heterozygous differences. Remaining differences are used to calculate quality scores.
