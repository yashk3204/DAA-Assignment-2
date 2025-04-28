input_file = 'CA-HepTh.txt'
output_file = 'CA-HepTh-preprocessed.txt'

edges = []
node_set = set()

with open(input_file, 'r') as f:
    first_line = f.readline()
    for line in f:
        parts = line.strip().split()
        a, b = parts
        if a == b:
            continue
        edges.append((a, b))
        node_set.add(a)
        node_set.add(b)

node_list = sorted(list(node_set))
node_id_map = {node: idx for idx, node in enumerate(node_list)}

edge_set = set()
for a, b in edges:
    u, v = node_id_map[a], node_id_map[b]
    if u > v:
        u, v = v, u
    edge_set.add((u, v))

with open(output_file, 'w') as f:
    f.write(f"{len(node_list)} {len(edge_set)}\n")
    for u, v in sorted(edge_set):
        f.write(f"{u+1} {v+1}\n")
