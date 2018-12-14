from indra.util import read_unicode_csv, write_unicode_csv
from indra.databases import hgnc_client

updated_lines = []
for line in read_unicode_csv('curated_site_map.csv', delimiter=','):
    gene = line[0]
    if gene == 'Gene':
        up_id = 'UniprotId'
    else:
        hgnc_id = hgnc_client.get_hgnc_id(gene)
        up_id = hgnc_client.get_uniprot_id(hgnc_id)
    updated_line = [up_id] + line
    updated_lines.append(updated_line)

write_unicode_csv('curated_site_map_up.csv', updated_lines)

