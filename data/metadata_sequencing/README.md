To generate those files

I go from the all_samples.tsv files needed for the rapidrun pipeline 

I select only the sample_name & Run columns (in this order)
You need the column names!

Then, I export as .csv with ';' separator

As it is standardised to be a .table file with a ' ' separator, I run the following command to change the extension+change the ';' to ' '. 

# Ex with Eparses file
tr ';' ' ' <Teleo_Eparses_clean.csv > Teleo_Eparses_clean.table

