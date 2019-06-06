import json
from flask import Flask, request, abort, Response, jsonify
from flask_cors import CORS
from protmapper import ProtMapper


app = Flask(__name__)
CORS(app)

optional_bool_args = {
    'do_methionine_offset': True,
    'do_orthology_mapping': True,
    'do_isoform_mapping': True
}
pm = ProtMapper()


@app.route('/map_to_human_ref', methods=['GET', 'POST'])
def map_to_human_ref():
    required_args = ('prot_id', 'prot_ns', 'residue', 'position')

    # Require all required arguments
    for arg in required_args:
        if request.json.get(arg) is None:
            abort(Response('The required argument "%s" is missing.' % arg,
                           400))
    # Now set the required arguments and the optional ones with default values
    # as backup
    arg_values = {key: request.json.get(key) for key in required_args}
    for arg, default_value in optional_bool_args.items():
        value = request.json.get(arg, default_value)
        arg_values[arg] = value

    ms = pm.map_to_human_ref(**arg_values)
    return jsonify(ms.to_json())


@app.route('/map_sitelist_to_human_ref', methods=['GET', 'POST'])
def map_sitelist_to_human_ref():
    site_list = request.json.get('site_list')
    if site_list is None:
            abort(Response('The required site_list argument is missing.', 400))

    for site in site_list:
        if len(site) != 4:
            abort(Response('Site list entries need to have exactly 4 elements.',
                           400))

    arg_values = {'site_list': site_list}
    for arg, default_value in optional_bool_args.items():
        value = request.json.get(arg, default_value)
        arg_values[arg] = value

    ms_list = pm.map_sitelist_to_human_ref(**arg_values)
    return jsonify([ms.to_json() for ms in ms_list])
