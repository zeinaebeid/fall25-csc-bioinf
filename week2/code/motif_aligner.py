from typing import List, Tuple, Optional, Dict

class MotifAligner:
    def align(
        self,
        sample_ids: List[str],
        encoded_vntrs: List[str],
        vid: Optional[str] = None,
        score_matrix: Optional[Dict] = None,
        output_dir: str = "./",
        tool: str = "center_star",
    ) -> Tuple[List[str], List[str]]:
        if tool == "center_star":
            return self._align_motifs_with_center_star(sample_ids, encoded_vntrs)
        raise ValueError(f"Unsupported alignment tool: {tool}")

    @staticmethod
    def _align_motifs_with_center_star(
        sample_ids: List[str],
        encoded_vntrs: List[str],
    ) -> Tuple[List[str], List[str]]:
        """
        Simple, deterministic 'center-star' placeholder:
        - Use the first sequence as the initial center.
        - For each next sequence, if it is longer than the current max length,
          pad ALL previously stored aligned sequences to the new max.
        - Append the new sequence padded to the current max.
        This guarantees all aligned strings end up with identical lengths.
        """
        if not encoded_vntrs:
            return sample_ids, encoded_vntrs

        aligned_ids: List[str] = []
        aligned_seqs: List[str] = []
        current_max = 0

        for sid, seq in zip(sample_ids, encoded_vntrs):
            if len(seq) > current_max:
                # New longest: pad everything we already added to this new length
                current_max = len(seq)
                for i in range(len(aligned_seqs)):
                    aligned_seqs[i] = aligned_seqs[i] + "-" * (current_max - len(aligned_seqs[i]))
            # Pad this sequence up to current_max and append
            padded = seq + "-" * (current_max - len(seq))
            aligned_ids.append(sid)
            aligned_seqs.append(padded)

        return aligned_ids, aligned_seqs

    @staticmethod
    def load_aligned_trs(aln_output: str) -> Tuple[List[str], List[str]]:
        aligned_trs: List[str] = []
        aligned_sample_ids: List[str] = []
        try:
            with open(aln_output, "r") as handle:
                for line in handle:
                    if line.startswith(">"):
                        aligned_sample_ids.append(line[1:].strip())
                    else:
                        aligned_trs.append(line.strip())
        except Exception:
            pass
        return aligned_sample_ids, aligned_trs
