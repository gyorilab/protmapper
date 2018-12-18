import json
from flask import Flask, request, abort, Response
from flask_cors import CORS
from sitemapper import SiteMapper


app = Flask(__name__)
CORS(app)


@app.route('/map/to_human_ref', methods=['GET', 'POST'])
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

    sm = SiteMapper()
    ms = sm.map_to_human_ref(**arg_values)
    return Response(json.dumps(ms.to_json()), mimetype='application/json')


if __name__ == '__main__':
    app.run()
