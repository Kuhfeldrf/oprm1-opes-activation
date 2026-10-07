"""Runbook §6.5: assert canonical P35372 numbering, author chain IDs, gaps,
SEQADV (engineered mutations) and non-receptor content for every candidate.

Usage: python cv/assert_numbering.py structures/raw
Exits non-zero on any numbering mismatch."""
import sys
from pathlib import Path
import gemmi

TARGETS = {166: "ASP", 167: "ARG", 168: "TYR", 281: "THR",
           334: "ASN", 335: "PRO", 336: "VAL", 337: "LEU", 338: "TYR", 339: "ALA"}
EXPECTED_CHAINS = {"10TM": {"R", "F"}, "8F7Q": {"R", "M"}, "8Y72": {"R"},
                   "9PPQ": {"R"}, "9PXU": {"R"}, "9MQH": {"A"}, "9MQI": {"A"},
                   "9MQJ": {"A"}}
IDS = ["9MQH", "9MQJ", "9MQI", "9PXU", "10TM", "9PPQ", "8F7Q", "8Y72"]

raw = Path(sys.argv[1] if len(sys.argv) > 1 else "structures/raw")
failed = False
for pid in IDS:
    doc = gemmi.cif.read(str(raw / f"{pid}.cif"))
    blk = doc.sole_block()
    st = gemmi.make_structure_from_block(blk)
    st.setup_entities()
    title = blk.find_value("_struct.title")
    res = blk.find_value("_em_3d_reconstruction.resolution") or blk.find_value("_refine.ls_d_res_high")
    print(f"== {pid}  res={res}  title={title}")
    # DBREF
    t = blk.find("_struct_ref_seq.", ["pdbx_strand_id", "pdbx_db_accession",
                 "seq_align_beg", "seq_align_end", "db_align_beg", "db_align_end",
                 "pdbx_auth_seq_align_beg", "pdbx_auth_seq_align_end"])
    receptor_chains = set()
    for r in t:
        print("   DBREF", " ".join(r))
        if r[1] == "P35372":
            receptor_chains.add(r[0])
    # SEQADV
    t = blk.find("_struct_ref_seq_dif.", ["pdbx_pdb_strand_id", "mon_id",
                 "pdbx_auth_seq_num", "db_mon_id", "pdbx_seq_db_seq_num", "details"])
    for r in t:
        if r[1] != "?" and r[3] != "?" and r[1] != r[3]:
            print("   SEQADV", " ".join(r))
    # entities
    for ent in st.entities:
        print(f"   entity {ent.name:>3} {ent.entity_type.name:<11} subchains={list(ent.subchains)}")
    for ch in st[0]:
        nums = {r.seqid.num: r.name for r in ch if r.het_flag == "A"}
        hits = {n: nums.get(n) for n in TARGETS if n in nums}
        if ch.name not in receptor_chains or not nums:
            hets = sorted({r.name for r in ch if r.het_flag == "H" and r.name != "HOH"})
            npoly = sum(1 for r in ch if r.het_flag == "A")
            print(f"   partner auth chain {ch.name}: {npoly} polymer residues, hetero {hets or 'none'}")
            continue
        bad = {n: v for n, v in hits.items() if v != TARGETS[n]}
        missing = [n for n in TARGETS if n not in nums]
        ok = not bad and not missing
        failed |= bool(bad) or ch.name not in EXPECTED_CHAINS[pid]
        print(f"   receptor auth chain {ch.name}: {'OK' if ok else f'MISMATCH {bad} missing {missing}'}"
              f"  span {min(nums)}-{max(nums)}  chain_expected={ch.name in EXPECTED_CHAINS[pid]}")
        prot = sorted(n for n in nums if 1 <= n <= 400)
        gaps = [n for n in range(prot[0], prot[-1] + 1) if n not in nums]
        print("   gaps (receptor):", gaps or "none")
        hets = sorted({r.name for r in ch if r.het_flag == "H" and r.name != "HOH"})
        print("   hetero in receptor chain:", hets or "none")
sys.exit(1 if failed else 0)
