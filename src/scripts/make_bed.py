import re
input_file = snakemake.input[0]
output_file = snakemake.output[0]

with open(input_file, "r") as input_file, open(output_file, "w") as out_file:
	header = next(input_file)
	header = header.strip().split("\t")
	count = 0
	for line in input_file:
		line = line.strip().split("\t")
		full_line = dict(zip(header, line))
		if float(full_line["log2fc"]) > 0:
			out_file.write(full_line["chromosome"] + "\t" + full_line["start"] + "\t" +
				full_line["end"] + "\t.\t" + full_line["log2fc"] + "\t" + full_line["strand"] + "\t" +
				full_line["peak_counts_clip"] + "\t" + full_line["pvalue"] + "\t" +
				full_line["total_clip_counts"] + "\n")
