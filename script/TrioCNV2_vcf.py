import sys

def process_trioCNV2(input_txt, output_vcf):
    """Process TrioCNV2 TXT directly to filtered, annotated VCF"""
    with open(input_txt, 'r') as f_in, open(output_vcf, 'w') as f_out:
        # Write VCF header
        f_out.write("##fileformat=VCFv4.3\n")
        f_out.write("##source=trioCNV2_processed\n")
        f_out.write('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\n')
        f_out.write('##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant">\n')
        f_out.write('##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of the variant">\n')
        f_out.write('##INFO=<ID=INHERITANCE,Number=1,Type=String,Description="Inheritance pattern: AD|AR|XL|YL|DeNovo|UPD|Unknown">\n')
        f_out.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        f_out.write('##FORMAT=<ID=CN,Number=1,Type=Integer,Description="Copy number">\n')
        f_out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG003NA24149Father\tHG004NA24143Mother\tHG002NA24385Son\n")

        for line in f_in:
            if line.startswith("#CHROM"):
                continue  # Skip header line
            
            # Skip empty lines
            if not line.strip():
                continue
                
            fields = line.strip().split()
            if len(fields) < 8:  # Ensure we have all expected columns
                continue
                
            chrom, start, end, father, mother, son, qual, evidence = fields[:8]

            # Determine SV type and ALT
            svtype = None
            alt = "."
            if "DELETION" in [father, mother, son]:
                svtype = "DEL"
                alt = "<DEL>"
            elif "DUPLICATION" in [father, mother, son]:
                svtype = "DUP"
                alt = "<DUP>"
            else:
                continue  # Skip non-SV lines

            # Calculate SV length
            svlen = int(end) - int(start)
            if svtype == "DEL":
                svlen = -svlen

            # Skip if SV length < 50 bp
            if abs(svlen) < 50:
                continue

            # Genotype conversion
            def get_gt(call):
                if call == "REFERENCE":
                    return "0/0:2"
                elif call == "DELETION":
                    return "0/1:1"  # Heterozygous deletion
                elif call == "DUPLICATION":
                    return "1/1:3"  # Homozygous duplication
                return "./.:."

            gt_father = get_gt(father)
            gt_mother = get_gt(mother)
            gt_son = get_gt(son)

            # Skip if both parents have variant but child doesn't
            if (gt_father != "0/0:2" and gt_mother != "0/0:2") and gt_son == "0/0:2":
                continue

            # Determine inheritance pattern
            def assign_inheritance(gt_f, gt_m, gt_c, chr):
                gt_f = gt_f.split(":")[0]
                gt_m = gt_m.split(":")[0]
                gt_c = gt_c.split(":")[0]

                # De Novo
                if (gt_c in ["0/1", "1/1"]) and (gt_f == "0/0" and gt_m == "0/0"):
                    return "DeNovo"
                
                # Autosomal Dominant
                if gt_c in ["0/1", "1/1"] and (gt_f != "0/0" or gt_m != "0/0"):
                    return "AD"
                
                # Autosomal Recessive
                if gt_c == "1/1" and gt_f == "0/1" and gt_m == "0/1":
                    return "AR"
                
                # UPD
                if gt_c == "1/1" and ((gt_f == "1/1" and gt_m == "0/0") or (gt_f == "0/0" and gt_m == "1/1")):
                    return "UPD"
                
                # X/Y-linked
                if chr == "chrX":
                    if gt_c == "1" and gt_m == "0/1":
                        return "XL"
                elif chr == "chrY":
                    if gt_f == "1" and gt_c == "1":
                        return "YL"
                
                return "Unknown"

            inheritance = assign_inheritance(gt_father, gt_mother, gt_son, chrom)

            # Write VCF line
            vcf_line = (
                f"{chrom}\t{start}\t.\tN\t{alt}\t{qual}\tPASS\t"
                f"SVTYPE={svtype};END={end};SVLEN={svlen};INHERITANCE={inheritance}\tGT:CN\t"
                f"{gt_father}\t{gt_mother}\t{gt_son}\n"
            )
            f_out.write(vcf_line)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python trioCNV2_processor.py <input.txt> <output.vcf>")
        sys.exit(1)
    
    input_txt = sys.argv[1]
    output_vcf = sys.argv[2]
    
    print(f"Processing {input_txt} directly to {output_vcf}...")
    process_trioCNV2(input_txt, output_vcf)
    print(f"Processing complete. Final output saved to {output_vcf}")