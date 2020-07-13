# coding=utf-8
import os
import share
import ntpath
import MINTplates


# getExpMatrixFile(blast_out_dir,out_dir)
# print("Matrix file created completely")

# tRNA_tRF_matrix: Sample/tRNA vs tRF_Type
# Sample_tRF_matrix: Sample vs tRF_Type
# tRNA_profile:Sample/tRNA/profile_type/profile_str
# ProfileType:start_pos, end_pos, uniq_profile, multiple_profile, total_profile

def getTrfReportFile(proj_name, dir, tRNA_dic, sample_dic, out_dir, offset=60):
    mean_trf_sample_dic = {}
    uniq_trf_sample_dic = {}
    trf_type_dic = {}
    mean_sample_trna_trftype_dic = {}
    uniq_sample_trna_trftype_dic = {}
    mean_sample_trftype_dic = {}
    uniq_sample_trftype_dic = {}
    mean_sample_trna_proftype_prof_dic = {}
    uniq_sample_trna_proftype_prof_dic = {}

    statistic_dic = {}

    sample_ls = []
    rna_ls = list(tRNA_dic.keys())
    trf_type_ls = []
    profile_type_ls = ["start_pos", "end_pos", "total"]

    mean_trna_trf_sample_matrix_file = out_dir + "/trf_sample_matrix.tsv"  # trf vs sample matrix
    MEAN_RNA_TRF_SAMPLE_MATRIX = open(mean_trna_trf_sample_matrix_file, 'w')


    mean_trna_trf_matrix_file = out_dir + "/trna_trftype_matrix.tsv"  # sample:trna vs trf type matrix
    MEAN_RNA_TRF_MATRIX = open(mean_trna_trf_matrix_file, 'w')  # For TRF type

    mean_sample_trf_matrix_file = out_dir + "/sample_trftype_matrix.tsv"  # Sample vs trf types (12) matrix
    MEAN_SAMPLE_TRF_MATRIX = open(mean_sample_trf_matrix_file, 'w')

    mean_profile_file = out_dir + "/profiles.tsv"  # Profiles tab files (including starpos, endpos, total)
    MEAN_PROFILE = open(mean_profile_file, 'w')

    combined_cleavage_sites = out_dir + "/cleavage_sites.tsv"  # Profiles tab files (including starpos, endpos, total)
    MEAN_CLEAVAGE = open(combined_cleavage_sites, 'w')

    trna_sample_readcount_matrix_file = out_dir + "/trna_sample_readcount_matrix.tsv"  # trna vs sample: read_count matrix
    SAMPLE_ReadCount_MATRIX = open(trna_sample_readcount_matrix_file, 'w')

    trna_sample_pileup_matrix_file = out_dir + "/trna_sample_pileup_matrix.tsv"  # trna vs sample: pileup max matrix
    SAMPLE_PileUP_MATRIX = open(trna_sample_pileup_matrix_file, 'w')

    hit_tab_ls = share.getExtFileList(dir, "hit.tab")

    # breif mapping informations for each hit
    sample_hit_proifle_dic = {}
    # For all sample
    for tab in hit_tab_ls:
        filenanme = ntpath.basename(tab)
        sample_id = os.path.basename(tab).replace("_" + proj_name + "_hit.tab", "")
        if sample_id in sample_dic:
            print("Preparing report for " + tab)
            statistic_dic[sample_id] = {
                "A": 0,
                "B": 0,
                "C": 0,
                "D": 0,
                "E": 0,
                "F": 0,
                "G": 0,
                "H": 0,
                "I": 0,
                "total":0,
                "u3_cl_ratio": 0,
                "intro_cl_ratio": 0,
                "u5_cl_ratio": 0,
                "cca_add_ratio": 0,
                "trna_mapped_read": 0
            }
            # sample_id = filenanme.replace("_"+proj_name+"_hit.tab","")
            if sample_id not in sample_ls:
                sample_ls.append(sample_id)
            if sample_id not in mean_sample_trna_trftype_dic:
                mean_sample_trna_trftype_dic[sample_id] = {}
            # if sample_id not in uniq_sample_trna_trftype_dic:
            #     uniq_sample_trna_trftype_dic[sample_id] = {}
            if sample_id not in mean_sample_trftype_dic:
                mean_sample_trftype_dic[sample_id] = {}
            # if sample_id not in uniq_sample_trftype_dic:
            #     uniq_sample_trftype_dic[sample_id] = {}
            if sample_id not in mean_sample_trna_proftype_prof_dic:
                mean_sample_trna_proftype_prof_dic[sample_id] = {}
            # if sample_id not in uniq_sample_trna_proftype_prof_dic:
            #     uniq_sample_trna_proftype_prof_dic[sample_id] = {}
            TAB = open(tab, 'r')

            # profile_dic
            # For each TRF
            trna_brief_infor_dic = {}
            for tab_line in TAB:
                if not tab_line.startswith("#"):
                    contents = tab_line.strip().split("\t")
                    rna_id = contents[1]
                    rna_family = contents[0]

                    if rna_id not in rna_ls:
                        rna_ls.append(rna_id)
                    # unimulti_total = float(contents[3])
                    if len(contents) < 29:
                        print(contents)
                    mean_exp = float(contents[28])
                    read_type = contents[32]
                    for a in read_type:
                        if a in statistic_dic[sample_id]:
                            statistic_dic[sample_id][a] += mean_exp
                            statistic_dic[sample_id]["total"] += mean_exp
                    statistic_dic[sample_id]["trna_mapped_read"] += mean_exp
                    hit_trna_number = int(contents[29])  # Unique/Multiple 1 or above 1
                    trf_type = contents[30]  # Nine Types of tRFs
                    if trf_type not in trf_type_ls:
                        trf_type_ls.append(trf_type)

                    brief_infor = contents[31]
                    sseq = brief_infor.split(",")[5]

                    tRF_id = MINTplates.encode_sequence(sseq, "")
                    if tRF_id not in trf_type_dic:
                        trf_type_dic[tRF_id] = {'seq': sseq, 'trf_type': trf_type, 'hit_trna_number': hit_trna_number,
                                                'tRNA_families': [rna_family], 'tRNA_IDs': [rna_id]}
                    else:
                        obj = trf_type_dic[tRF_id]
                        if rna_family not in obj['tRNA_families']:
                            obj['tRNA_families'].append(rna_family)
                        if rna_id not in obj['tRNA_IDs']:
                            obj['tRNA_IDs'].append(rna_id)

                    if tRF_id not in mean_trf_sample_dic:
                        mean_trf_sample_dic[tRF_id] = {}
                    # if tRF_id not in uniq_trf_sample_dic:
                    #     uniq_trf_sample_dic[tRF_id] = {}
                    # if sample_id not in uniq_trf_sample_dic[tRF_id]:
                    #     uniq_trf_sample_dic[tRF_id][sample_id]={'value': 0}
                    if sample_id not in mean_trf_sample_dic[tRF_id]:
                        mean_trf_sample_dic[tRF_id][sample_id] = {'value': 0}
                    mean_trf_sample_dic[tRF_id][sample_id]['value'] += mean_exp

                    if trf_type not in mean_sample_trftype_dic[sample_id]:
                        mean_sample_trftype_dic[sample_id][trf_type] = 0

                    if rna_id not in mean_sample_trna_trftype_dic[sample_id]:
                        mean_sample_trna_trftype_dic[sample_id][rna_id] = {}

                    if trf_type not in mean_sample_trna_trftype_dic[sample_id][rna_id]:
                        mean_sample_trna_trftype_dic[sample_id][rna_id][trf_type] = 0

                    mean_sample_trftype_dic[sample_id][trf_type] += mean_exp
                    # uniq_sample_trftype_dic[sample_id][trf_type]+=unimulti_total
                    mean_sample_trna_trftype_dic[sample_id][rna_id][trf_type] += mean_exp
                    # uniq_sample_trna_trftype_dic[sample_id][rna_id][trf_type]+=unimulti_total

                    # Profile
                    if rna_id not in trna_brief_infor_dic:
                        trna_brief_infor_dic[rna_id] = []
                    trna_brief_infor_dic[rna_id].append(brief_infor)
            # Calculate cleavage ratio for each samples
            obj = statistic_dic[sample_id]
            if obj["H"] + obj["G"] + obj["I"] > 0:
                obj["u3_cl_ratio"] = round((obj["H"] + obj["I"]) / (obj["H"] + obj["G"] + obj["I"]), 3)
            if obj["E"] + obj["F"] > 0:
                obj["intro_cl_ratio"] = round(obj["F"] / (obj["E"] + obj["F"]), 3)
            if obj["C"] + obj["B"] > 0:
                obj["u5_cl_ratio"] = round((obj["C"]) / (obj["C"] + obj["B"]), 3)
            if obj["I"] + obj["H"] > 0:
                obj["cca_add_ratio"] = round((obj["I"]) / (obj["I"] + obj["H"]), 3)

            # Get profile and cleavage sites informataion
            profile_dic = {}
            for rna_id in trna_brief_infor_dic:
                if rna_id in tRNA_dic:
                    t = tRNA_dic[rna_id]
                    profile_dic[rna_id] = t.CreateTrfProfiles(trna_brief_infor_dic[rna_id], offset)
            if sample_id not in sample_hit_proifle_dic:
                sample_hit_proifle_dic[sample_id] = profile_dic

    # Print expression matrix of tRFs
    MEAN_RNA_TRF_SAMPLE_MATRIX.write("tRF_ID\ttRNA_Families\tRNA_IDs\tseq\tType\t" + "\t".join(sample_ls) + "\n")
    for trf_id in mean_trf_sample_dic:
        seq = trf_type_dic[trf_id]['seq']
        c_trf_type = trf_type_dic[trf_id]['trf_type']
        # hit_trna_number = trf_type_dic[trf_id]['hit_trna_number']
        tRNA_Families = ",".join(trf_type_dic[trf_id]['tRNA_families'])
        tRNA_IDs = ",".join(trf_type_dic[trf_id]['tRNA_IDs'])

        line = trf_id + "\t" + tRNA_Families + "\t" + tRNA_IDs + "\t" + seq + "\t" + c_trf_type
        for sample_id in sample_ls:
            value = 0
            if sample_id in mean_trf_sample_dic[trf_id]:
                value = round(mean_trf_sample_dic[trf_id][sample_id]["value"], 3)
            line += "\t" + str(value)
        MEAN_RNA_TRF_SAMPLE_MATRIX.write(line + "\n")
    MEAN_RNA_TRF_SAMPLE_MATRIX.close()

    # UNIQ_RNA_TRF_SAMPLE_MATRIX.write("tRF_ID\tRNA_family\tRNA_ID\ttRF_ID\tType\tHit_Type\t" + "\t".join(sample_ls) + "\n")
    # for trf_id in uniq_trf_sample_dic:
    #     c_trf_type = trf_type_dic[trf_id]['trf_type']
    #     c_hit_type = trf_type_dic[trf_id]['hit_type']
    #     tRNA_id = trf_id.split("#")[0]
    #     family_id = tRNA_id
    #     if tRNA_id in tRNA_id_family_dic:
    #         family_id = tRNA_id_family_dic[tRNA_id]
    #     line = trf_id+"\t"+family_id + "\t" + "\t".join(trf_id.split("#")) + "\t" + c_trf_type + "\t" + c_hit_type
    #     for sample_id in sample_ls:
    #         value = 0
    #         if sample_id in uniq_trf_sample_dic[trf_id]:
    #             value = uniq_trf_sample_dic[trf_id][sample_id]["value"]
    #         line += "\t" + str(value)
    #     UNIQ_RNA_TRF_SAMPLE_MATRIX.write(line + "\n")
    # UNIQ_RNA_TRF_SAMPLE_MATRIX.close()

    MEAN_RNA_TRF_MATRIX.write("#SampleID\tRNA_family\tRNA_ID\t" + "\t".join(trf_type_ls) + "\n")
    # UNIQ_RNA_TRF_MATRIX.write("#SampleID\tRNA_family\tRNA_ID\t" + "\t".join(trf_type_ls) + "\n")
    for sample in sample_ls:
        for rna_id in rna_ls:
            if rna_id in tRNA_dic:
                rna = tRNA_dic[rna_id]
                mean_line = sample + "\t" + rna.family + "\t" + rna_id
                for trf_type in trf_type_ls:
                    if sample in mean_sample_trna_trftype_dic \
                            and rna_id in mean_sample_trna_trftype_dic[sample] \
                            and trf_type in mean_sample_trna_trftype_dic[sample][rna_id]:
                        mean_line += "\t" + str(round(mean_sample_trna_trftype_dic[sample][rna_id][trf_type], 3))
                    else:
                        mean_line += "\t0"
                MEAN_RNA_TRF_MATRIX.write(mean_line + "\n")
    MEAN_RNA_TRF_MATRIX.close()

    MEAN_SAMPLE_TRF_MATRIX.write("#SampleID\t" + "\t".join(trf_type_ls) + "\n")
    for sample in sample_ls:
        mean_line = sample
        for trf_type in trf_type_ls:
            if sample in mean_sample_trftype_dic \
                    and trf_type in mean_sample_trftype_dic[sample]:
                mean_line += "\t" + str(round(mean_sample_trftype_dic[sample][trf_type], 3))
            else:
                mean_line += "\t0"
        MEAN_SAMPLE_TRF_MATRIX.write(mean_line + "\n")
    MEAN_SAMPLE_TRF_MATRIX.close()

    SAMPLE_ReadCount_MATRIX.write("#tRNA_family\ttRNA_ID\t" + "\t".join(sample_ls)+"\n")
    SAMPLE_PileUP_MATRIX.write("#tRNA_family\ttRNA_ID\t" + "\t".join(sample_ls)+"\n")
    for rna_id in rna_ls:
        if rna_id in tRNA_dic:
            rna = tRNA_dic[rna_id]
            family = rna.family
            rc_line = family + "\t" + rna_id
            ph_line = family + "\t" + rna_id
            for sample in sample_ls:
                total_reads_num = 0
                pileup_height = 0
                if sample in sample_hit_proifle_dic \
                        and rna_id in sample_hit_proifle_dic[sample]:
                    total_reads_num = sample_hit_proifle_dic[sample][rna_id]["total_reads_num"]
                    pileup_height = sample_hit_proifle_dic[sample][rna_id]["pileup_height"]
                rc_line += "\t" + str(round(total_reads_num,3))
                ph_line += "\t" + str(round(pileup_height,3))
            SAMPLE_ReadCount_MATRIX.write(rc_line + "\n")
            SAMPLE_PileUP_MATRIX.write(ph_line + "\n")
    SAMPLE_ReadCount_MATRIX.close()
    SAMPLE_PileUP_MATRIX.close()

    MEAN_PROFILE.write("#SampleID\ttRNA_family\ttRNA_ID\ttype\tprofile\n")
    MEAN_CLEAVAGE.write(
        "#SampleID\ttRNA_family\ttRNA_ID\tID\tPosition\tPType\tIntensity\tIntensity_5\tIntensity_3\tSNRatio\tSeq_5\tSeq_3\tfull_Seq\n")
    for sample in sample_ls:
        for rna_id in rna_ls:
            family = "Unknown"
            if rna_id in tRNA_dic:
                rna = tRNA_dic[rna_id]
                family = rna.family
            for prof_type in profile_type_ls:
                if sample in sample_hit_proifle_dic \
                        and rna_id in sample_hit_proifle_dic[sample] \
                        and prof_type in sample_hit_proifle_dic[sample][rna_id]:
                    mean_line = sample + "\t" + family + "\t" + rna_id + "\t" + prof_type + "\t" + ",".join(
                        share.stringfyNumList(sample_hit_proifle_dic[sample][rna_id][prof_type]))
                    MEAN_PROFILE.write(mean_line + "\n")

            if sample in sample_hit_proifle_dic \
                and rna_id in sample_hit_proifle_dic[sample] \
                and prof_type in sample_hit_proifle_dic[sample][rna_id]:
                site_dic = sample_hit_proifle_dic[sample][rna_id]["cleavage_site_dic"]
                for site_id in site_dic:
                    site_obj = site_dic[site_id]
                    cleavage_site_str = sample + "\t" + family + "\t" + rna_id + \
                                        "\t" + str(site_id) + \
                                        "\t" + str(site_obj["pos"]) + \
                                        "\t" + str(site_obj["ptype"]) + \
                                        "\t" + str(int(site_obj["int"])) + \
                                        "\t" + str(int(site_obj["int5"])) + \
                                        "\t" + str(int(site_obj["int3"])) + \
                                        "\t" + str(round(float(site_obj["sn_ratio"]), 3)) + \
                                        "\t" + site_obj["c_5_seq"] + \
                                        "\t" + site_obj["c_3_seq"] + \
                                        "\t" + site_obj["c_5_seq"] + site_obj["c_3_seq"]
                    MEAN_CLEAVAGE.write(cleavage_site_str + "\n")
    print("Output 1: " + mean_trna_trf_sample_matrix_file)
    print("Output 2: " + mean_trna_trf_matrix_file)
    print("Output 3: " + mean_sample_trf_matrix_file)
    print("Output 4: " + mean_profile_file)
    print("Output 5: " + combined_cleavage_sites)
    print("Output 6: " + trna_sample_readcount_matrix_file)
    print("Output 7: " + trna_sample_pileup_matrix_file)

    MEAN_PROFILE.close()
    MEAN_CLEAVAGE.close()
    return statistic_dic
