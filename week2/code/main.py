from typing import List, Optional, Tuple
from decomposer import Decomposer
from motif_encoder import MotifEncoder
from motif_aligner import MotifAligner

class TandemRepeatVizWorker:
    def __init__(self):
        self.decomposer = Decomposer()
        self.motif_encoder = MotifEncoder()
        self.motif_aligner = MotifAligner()

    def generate_trplot(
        self,
        tr_id: str,
        sample_ids: List[str],
        tr_sequences: List[str],
        motifs: List[str],
        skip_alignment: bool = False,
        rearrangement_method: str = "motif_count",
        sample_order_file: Optional[str] = None,   # ✅ Optional
        output_dir: str = "./",
        verbose: bool = True,
    ) -> Tuple[List[str], List[str], List[List[str]]]:
        """
        Main workflow: decomposition → encoding → alignment → rearrangement
        """
        if len(sample_ids) != len(tr_sequences):
            raise ValueError("The number of sample IDs and sequences are different.")

        if verbose:
            print(f"[worker] ID: {tr_id}")
            print(f"[worker] Motifs: {motifs}")
            print(f"[worker] Loaded {len(tr_sequences)} TR sequences")
            print("[worker] Decomposing…")

        # 1) Decomposition
        decomposed_trs: List[List[str]] = []
        for i, tr_sequence in enumerate(tr_sequences):
            if verbose:
                print(f"[worker]  progress {i+1}/{len(tr_sequences)}")
            decomposed_trs.append(self.decomposer.decompose(tr_sequence, motifs))

        # 2) Refinement
        decomposed_trs = self.decomposer.refine(decomposed_trs, verbose=False)

        # 3) Encoding
        if verbose:
            print("[worker] Encoding")
        encoded_trs = self.motif_encoder.encode(
            decomposed_trs,
            motif_map_file=f"{output_dir}/{tr_id}_motif_map.txt",
            auto=True,
        )

        # 4) Alignment
        if skip_alignment:
            if verbose:
                print("[worker] Skipping alignment (pad)")
            max_len = max(len(seq) for seq in encoded_trs) if encoded_trs else 0
            aligned_trs = [seq + "-" * (max_len - len(seq)) for seq in encoded_trs]
            sorted_sample_ids = sample_ids
        else:
            if verbose:
                print("[worker] Alignment (center-star)")
            sorted_sample_ids, aligned_trs = self.motif_aligner.align(
                sample_ids, encoded_trs, vid=tr_id, output_dir=output_dir, tool="center_star"
            )

        if verbose:
            print("[worker] Done")

        return sorted_sample_ids, aligned_trs, decomposed_trs
