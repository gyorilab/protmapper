import json
from flask import Flask, request, abort, Response, jsonify
from flask_cors import CORS
from protmapper import ProtMapper


app = Flask(__name__)
CORS(app)


@app.route('/map_to_human_ref', methods=['GET', 'POST'])
def map_to_human_ref():
    args = ('prot_id', 'prot_ns', 'residue', 'position',
            'do_methionine_offset', 'do_orthology_mapping',
            'do_isoform_mapping')
    required_args = args[:4]
    for arg in required_args:
        if request.json.get(arg) is None:
            abort(Response('The required argument "%s" is missing.' % arg,
                           400))

    arg_values = {key: request.json.get(key) for key in args}

    pm = ProtMapper()
    ms = pm.map_to_human_ref(**arg_values)
    return jsonify(ms.to_json())


@app.route('/map_sitelist_to_human_ref', methods=['GET', 'POST'])
def map_sitelist_to_human_ref():
    site_list = request.json.get('site_list')
    if site_list is None:
            abort(Response('The required sute_list argument is missing.' %
                           arg, 400))

    for site in site_list:
        if len(site) != 4:
            abort(Response('Site list entries need to have exactly '
                           '4 elements.', 400))

    opt_arg_values = {key: request.json.get(key) for key in
                      ('do_methionine_offset', 'do_orthology_mapping',
                       'do_isoform_mapping')}
    pm = ProtMapper()
    ms_list = pm.map_sitelist_to_human_ref(site_list=site_list,
                                           **opt_arg_values)
    return jsonify([ms.to_json() for ms in ms_list])

if __name__ == '__main__':
    app.run()
