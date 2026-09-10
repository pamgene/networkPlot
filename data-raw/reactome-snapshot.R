# Bundled Reactome reference snapshot in inst/extdata/
#
#   ReactomePathways.txt          (id \t name \t organism, no header)
#   ReactomePathwaysRelation.txt  (parent \t child, no header)
#
# These are shipped so reactome_refs() / enrich_network() work out of the
# box. To refresh when Reactome cuts a new release (~quarterly):
#
#   download_reactome_refs("networkPlot/inst/extdata")
#
# which fetches https://reactome.org/download/current/ReactomePathways.txt
# and .../ReactomePathwaysRelation.txt. reactome_refs() filters to
# "Homo sapiens" at load time.
#
# This snapshot was taken from Reactome's current release (copied from the
# Network_generation repo's data/ folder, dated 2025-11-24).
