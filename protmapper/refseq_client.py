import logging
from protmapper.resources import resource_manager
from protmapper.uniprot_client import load_fasta_sequences


logger = logging.getLogger(__name__)


class RefseqMapper(object):
    def __init__(self):
        self.initialized_ids = False
        self.initialized_seq = False

    def initialize_ids(self):
        self._refseq_uniprot = _build_refseq_entries()
        self.initialized_ids = True

    def initialize_seq(self):
        self._sequences = _build_refseq_sequences()
        self.initialized_seq = True

    @property
    def refseq_uniprot(self):
        if not self.initialized_ids:
            self.initialize_ids()
        return self._refseq_uniprot

    @property
    def sequences(self):
        if not self.initialized_seq:
            self.initialize_seq()
        return self._sequences


rm = RefseqMapper()


def _build_refseq_sequences():
    seq_file = resource_manager.get_create_resource_file('refseq_seq',
                                                         cached=True)
    logger.info("Loading RefSeq protein sequences...")
    seq = load_fasta_sequences(seq_file, id_delimiter=' ', id_index=0)
    return seq
