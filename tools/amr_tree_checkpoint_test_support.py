"""Locate tree checkpoint records for targeted corruption in restart checks."""


def checkpoint_records(data):
    offset = 0

    def line():
        nonlocal offset
        end = data.index(b'\n', offset)
        value = data[offset:end]
        offset = end + 1
        return value

    tag = line()
    if tag in (b'PELEF_AMR_TREE_INTERNAL_ENERGY_V1',
               b'PELEF_AMR_TREE_INTERNAL_ENERGY_V2'):
        line()  # Energy correction ledger precedes the checkpoint schema.
        tag = line()
    history = tag == b'PELEF_AMR_TREE_TIMESTEP_HISTORY_V1'
    if history:
        tag = line()
    # All six schemas share the context, metadata and boundary-ledger prefix.
    # Reaction sources and heat release extend records after that ledger.
    assert tag in tuple(f'PELEF_AMR_TREE_3D_CHECKPOINT_V{version}'.encode()
                        for version in range(1, 7)), tag
    size = int(line())
    assert size > 0 and offset + size < len(data)
    offset += size
    records = data[offset:].splitlines(keepends=True)
    assert len(records[0].split()) == 5
    assert len(records[1].split()) == 6
    return data[:offset], records, 2 + int(history), tag
