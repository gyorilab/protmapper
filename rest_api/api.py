from flask import Flask, request, abort, Response
from flask_cors import CORS
from sitemapper import SiteMapper


app = Flask(__name__)
CORS(app)


@app.route('/map/to_human_ref', methods=['GET', 'POST'])
def map_to_human_ref()
    prot_id = request.json.get('prot_id')
    prot_ns = request.json.get('prot_ns')
    residue = request.json.get('residue')
    position = request.json.get('position')
    do_methionine_offset = request.json.get('do_methionine_offset')
    do_orthology_mapping = request.json.get('do_orthology_mapping')
    do_isoform_mapping = request.json.get('do_isoform_mapping')

    sm = SiteMapper()
    ms = sm.map_to_human_ref(prot_id, prot_ns, residue, position,
                             do_methionine_offset=do_methionine_offset,
                             do_orthology_mapping=do_orthology_mapping,
                             do_isoform_mapping=do_isoform_mapping)
    return Response()
