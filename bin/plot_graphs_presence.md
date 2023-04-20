```
API = "https://www.ebi.ac.uk/metagenomics/api/v1"

def fetch(data, url, next_mode=True, logging_flag=False, return_field="attributes", params=None):
    logging.getLogger().setLevel(level=logging.DEBUG if logging_flag else logging.INFO) 
    logging.debug(f"fetching: {url}")
    response = requests.get(url, params=params)

    response_data = response.json()
    next_url = response_data.get("links", {}).get("next")

    data.extend([entry.get(return_field) for entry in response_data.get("data")])
    logging.debug(str(len(data)))
    if next_mode:
        if next_url:
            fetch(data, next_url, params, logging_flag, return_field)


def create_dot(name, presented):
    dot = graphviz.Digraph(name, comment=pathways_schema[name]) 
    edges = graph[0].edges
    max_weight = 0
    for edge, count in zip(edges, range(len(edges))):
        from_node = edge[0]
        to_node = edge[1]
        number = edge[2]
        dot.node(str(from_node))
        dot.node(str(to_node))

        label = edges._adjdict[from_node][to_node][number]['label']
        weight = edges._adjdict[from_node][to_node][number]['weight']
        if 1/weight > max_weight:
            max_weight = int(1/weight)
        if weight == 1 or weight == 0 :
            weight_str = str(weight)
        else:
            weight_str = '1/' + str(int(1/weight))
        color = 'red' if label in presented else 'black'
        dot.edge(str(from_node), str(to_node), label=label + ' \n [' + weight_str + ']', color=color)
    return dot
    
def generate_graph(tax, counts, nodes, source, target, value, edges):
    for line in tax:
        groups = line.split(':')
        if len(groups) == 1:
            name = levels[0] + '_' + groups[0]
            if name not in edges:
                source.append(0)
                target.append(nodes[name])
                value.append(counts[name])
                edges.append(name)
        else:                 
            for i in range(1, len(groups)):
                name = ';'.join([levels[j] + '_' + groups[j] for j in range(i+1)])
                if name not in edges:
                    parent_name = ';'.join([levels[j] + '_' + groups[j] for j in range(i)])
                    source.append(nodes[parent_name])
                    target.append(nodes[name])
                    value.append(counts[name])
                    edges.append(name)
    return source, target, value, edges


def create_nodes(nodes, tax, counts):
    counts['Life'] += sum(tax.values())
    number = len(nodes)
    for line in tax:
        groups = line.split(':')
        for i in range(len(groups)):
            name = ';'.join([levels[j] + '_' + groups[j] for j in range(i+1)])
            if name not in nodes:
                nodes.setdefault(name, number)
                number += 1
                counts.setdefault(name, 0)
        for i in range(len(groups)):
            name = ';'.join([levels[j] + '_' + groups[j] for j in range(i+1)])
            counts[name] += tax[line]
    return nodes, counts
    

pathways_schema = {}
with open('kegg/all_pathways.txt', 'r') as pathways_file:
    for line in pathways_file:
        fields = line.strip().split(':')
        pathways_schema[fields[0]] = fields[1]

downloads = []
fetch(data=downloads, 
      url=f"https://www.ebi.ac.uk/metagenomics/api/v1/analyses/{mgya}/downloads", 
      logging_flag=False, return_field='links')

for item in downloads:
    if item.get('self'):
        #if 'ko.tsv' in item.get('self'):
        if 'kegg_pathways.csv' in item.get('self'):
            ko_file = item.get('self')         
response = requests.get(ko_file)
data = response.text.split('\n')

value = randint(1, len(data)-1)
pathway = data[value].split('","')
presented_ko = pathway[4:-1][0].split(',')
name = pathway[0][1:]
print(f'Pathway: {name}, completeness: {pathway[1]}')
print(f'Presented KOs: {presented_ko}')

graph_file = open("kegg/graphs.pkl", 'rb')
graphs = pickle.load(graph_file)
graph = graphs[name]

dot = create_dot(name, presented=presented_ko)
dot.render(directory='kegg', filename=name, format='png')
Image(f'kegg/{name}.png')
```