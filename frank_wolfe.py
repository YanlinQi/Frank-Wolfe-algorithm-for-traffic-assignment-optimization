import networkx
import pandas as pd
# import osmnx as ox
import numpy as np


class Stack:
    def __init__(self):
        self.stack = []
        self.num_items = 0
        self.setc = set([])

    def get_next(self):
        next_itm = self.stack.pop()
        self.setc.remove(next_itm[1])
        self.num_items -= 1
        return next_itm

    def add_next(self, new_itm):
        self.setc.add(new_itm[1])
        self.stack.append(new_itm)
        self.num_items += 1

    def exists_in(self, itm):
        return itm in self.setc


def label_correcting(g, origin, destination, stack_obj, upper_init=np.inf):
    """
    find the shortest path of a directional graph using label_correcting shortest path algorithm
    :param g: graph object
    :param origin: origin node name on the graph
    :param destination: destination node name on the graph
    :param stack_obj: stack data structure for prioritizing search type
    :param upper_init:  initial maximum distance limit from origin to destination
    :return: the shortest path formed by node list;
             the shortest accumulated distance of the path;
             the node set that has been visited
    """
    stack_obj.add_next((0, origin))
    dist_shortest = {origin: {'dist': 0, 'parent': None}}
    upper = upper_init
    # keep track of all nodes visited
    nds_vsted = [origin]

    while stack_obj.num_items > 0:
        candi_node = stack_obj.get_next()
        for node in g.successors(candi_node[1]):
            # obtain distance from parent node to child node
            d_ij = g.get_edge_data(candi_node[1], node)['weight']
            # obtain distance to parent node
            d_i = candi_node[0]
            # distance to node j is the sum of distance to parent node i and the distance between node i and j
            d_j = d_i + d_ij
            # left hand side
            lhs = d_j

            if node in dist_shortest:
                # right hand side
                rhs = min(dist_shortest[node]['dist'], upper)
            else:
                rhs = min(np.inf, upper)

            if lhs < rhs:
                dist_shortest[node] = {'dist': d_j, 'parent': candi_node[1]}
                if (node != destination) and (not stack_obj.exists_in(node)):
                    stack_obj.add_next((d_j, node))
                    nds_vsted.append(node)

                if node == destination:
                    upper = d_j
    min_path_dist = dist_shortest[destination]['dist']
    r_path = [destination]
    next_node = destination
    flag = True
    # obtain the actual node_chain through back-tracking on dist_shortest
    while flag:
        next_node = dist_shortest[next_node]['parent']
        if next_node is None:
            flag = False
        else:
            r_path.append(next_node)
    # the actual path formed by the node chain
    path = list(reversed(r_path))

    return path, min_path_dist


class TrafficUE:
    def __init__(self):
        self.graph = networkx.DiGraph()
        self.od_demand = dict()

    def load_graph(self, g_file, col_names):
        """
        generate a bi-directional graph structure from the given .xlsx
        based on the orgin, destination, link-cost information
        :param g_file: .xlsx file including graph info
        :param col_names: [origin_column_name, destination_column_name, link-cost_column_name]
        :return: the generated graph structure of networkx.DiGraph()
        """
        net_df = pd.read_excel(g_file,
                               usecols=col_names,
                               dtype={col_names[0]: str, col_names[1]: str})

        net_info = net_df[col_names[0:-1]].to_records(index=False)
        self.graph.add_weighted_edges_from(net_info)
        for i, [o_i, d_i] in enumerate(net_info):
            self.graph[o_i][d_i]['x'] = 0.0
            self.graph[o_i][d_i]['y'] = 0.0
            self.graph[o_i][d_i]['t'] = self.graph[o_i][d_i]['y']-self.graph[o_i][d_i]['x']
            self.graph[o_i][d_i]['capacity'] = net_df[col_names[-1]][i]
            self.graph[o_i][d_i]['t0'] = self.graph[o_i][d_i]['weight']

        return self.graph

    def load_demand(self, od_file):
        with open(od_file, "r") as f:
            text = f.read()
            demand_info = list(map(lambda x: x.split(), text.splitlines()))
            num_origins = int(demand_info[0][0])
            od_ptr = list(map(lambda x: [x[0], int(x[1])], demand_info[1:num_origins + 1]))
            des_info = list(map(lambda x: [x[0], float(x[1])], demand_info[num_origins + 2:]))

            for i, [origin_i, des_start] in enumerate(od_ptr):
                if i < len(od_ptr) - 1:
                    des_end = od_ptr[i + 1][1] - 1
                else:
                    des_end = len(des_info)

                self.od_demand[origin_i] = dict()
                for des_j in range(des_start - 1, des_end):
                    self.od_demand[origin_i][des_info[j][0]] = des_info[des_j][1]
        return self.od_demand

    def update_link_cost(self, o, d, alpha=0.15, beta=4):
        self.graph[o][d]['weight'] = self.graph[o][d]['t0'] * \
                                             (1 + alpha * (pow((self.graph[o][d]['x'] /
                                                                (0.9 * self.graph[o][d]['capacity'])), beta)))

    def update_direction(self, o, d):
        self.graph[o][d]['d'] = self.graph[o][d]['y'] - self.graph[o][d]['x']

    def AoN(self, o, d):
        route_nodes, min_cost = label_correcting(self.graph, o, d, Stack())
        for num_node in range(len(route_nodes) - 1):
            self.graph[route_nodes[num_node]][route_nodes[num_node + 1]]['y'] = self.od_demand[o][d]


if __name__ == "__main__":
    traff = TrafficUE()
    traff.load_graph("anaheim.xlsx", ['start', 'end', 'Free Flow Travel Time (min)', 'Capacity'])
    traff.load_demand("fort2.txt")
    for o_i in traff.od_demand.keys():
        for d_i in traff.od_demand[o_i].keys():
            print('...')

