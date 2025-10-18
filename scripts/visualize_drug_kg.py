"""
Script: visualize_drug_kg.py
Purpose: Visualize the drug-target-cancer knowledge graph from Neo4j
Dependencies: neo4j, networkx, matplotlib
"""

from neo4j import GraphDatabase
import networkx as nx
import matplotlib.pyplot as plt

# ---------------------------
# Neo4j connection config
# ---------------------------
NEO4J_URI = "bolt://localhost:7687"
NEO4J_USER = "neo4j"
NEO4J_PASSWORD = "MyNewSecurePassword123"

driver = GraphDatabase.driver(NEO4J_URI, auth=(NEO4J_USER, NEO4J_PASSWORD))

# ---------------------------
# Fetch all KG relationships
# ---------------------------
def fetch_all_kg():
    query = """
    MATCH (n)
    OPTIONAL MATCH (n)-[r]->(m)
    RETURN n.name AS from_node, labels(n)[0] AS from_type,
           type(r) AS relation, m.name AS to_node, labels(m)[0] AS to_type
    """
    with driver.session() as session:
        result = session.run(query)
        return [r for r in result]


# ---------------------------
# Build NetworkX graph
# ---------------------------
def build_graph(data):
    G = nx.Graph()
    for record in data:
        from_node = record["from_node"]
        from_type = record["from_type"]
        to_node = record["to_node"]
        to_type = record["to_type"]

        # Add nodes with type
        if from_node:
            G.add_node(from_node, type=from_type.lower())
        if to_node:
            G.add_node(to_node, type=to_type.lower())

        # Add edge if there is a relation
        if record["relation"] and from_node and to_node:
            G.add_edge(from_node, to_node)

    return G


# ---------------------------
# Draw graph
# ---------------------------
def draw_graph(G):
    plt.figure(figsize=(12, 8))

    pos = nx.spring_layout(G, k=0.5, seed=42)
    
    # Node groups by type
    drugs = [n for n, attr in G.nodes(data=True) if attr["type"] == "drug"]
    targets = [n for n, attr in G.nodes(data=True) if attr["type"] == "target"]
    cancers = [n for n, attr in G.nodes(data=True) if attr["type"] == "cancer"]

    # Draw nodes
    nx.draw_networkx_nodes(G, pos, nodelist=drugs, node_color="skyblue", node_size=600, label="Drug")
    nx.draw_networkx_nodes(G, pos, nodelist=targets, node_color="lightgreen", node_size=500, label="Target")
    nx.draw_networkx_nodes(G, pos, nodelist=cancers, node_color="salmon", node_size=500, label="Cancer")

    # Draw edges
    nx.draw_networkx_edges(G, pos, width=1.0, alpha=0.7)

    # Draw labels
    nx.draw_networkx_labels(G, pos, font_size=9)

    plt.title("Drug-Target-Cancer Knowledge Graph")
    plt.legend(scatterpoints=1)
    plt.axis("off")
    plt.show()

# ---------------------------
# Main
# ---------------------------
if __name__ == "__main__":
    data = fetch_all_kg()
    if not data:
        print("No data found in the knowledge graph.")
    else:
        G = build_graph(data)
        draw_graph(G)
