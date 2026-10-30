# Bundled pathway data in inst/extdata/
#
#   reactome_v2_pathway_table.csv        pathway_id, pathway_name, uniprot: Reactome cut at a curated
#                                        depth per top-level category ("Reactome v2"), all human
#                                        UniProt members
#   reactome_v2_pathway_table_info.txt   Reactome release, download dates, depth table
#   reactome_v2_pathway_annotation.csv   per pathway: category, disease flags (is_in_disease,
#                                        disease_group, normal counterpart, excluded_by_default)
#   reactome_v2_normal_counterpart_members.csv  normal_pathway_id, uniprot: all-level members of
#                                        each normal counterpart, for filter_disease_pathways()
#   ukadb_0.6.0_kinases.csv              UKAdb 0.6.0 kinase list (Kinase_Name, Kinase_UniprotID, ...)
#
# Built in the Network_generation repo by
#   Rscript tools/build_reactome_pathway_table.R <reactome folder> [<Reactome release>]
# (Reactome flat files + Reactome ContentService for the disease flags) and copied here unchanged
# from Network_generation/data/. To update: rebuild there, then copy the five files into
# inst/extdata and re-run the tests.
