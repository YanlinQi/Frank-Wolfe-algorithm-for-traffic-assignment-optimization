"""
@Author: Yanlin Qi
Date: 5/10/2021
Implementation of Frank-Wolfe algorithm for traffic assignment optimizatoin
"""

import matplotlib.pyplot as plt
import networkx
import pandas as pd
import numpy as np
from scipy import integrate
import copy


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


class TrafficOptim:
    def __init__(self):
        self.graph = networkx.DiGraph()
        self.od_demand = dict()
        self.conv_df = pd.DataFrame()

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
        for i, [o_i, d_i, w_i] in enumerate(net_info):
            self.graph[o_i][d_i]['x'] = 0.0
            self.graph[o_i][d_i]['y'] = 0.0
            self.graph[o_i][d_i]['d'] = self.graph[o_i][d_i]['y'] - self.graph[o_i][d_i]['x']
            self.graph[o_i][d_i]['capacity'] = net_df[col_names[-1]][i]
            self.graph[o_i][d_i]['t0'] = self.graph[o_i][d_i]['weight']

        return self.graph

    def load_demand(self, od_file, scaling_coeff=1.0):
        with open(od_file, "r") as f:
            text = f.read()
            demand_info = list(map(lambda x: x.split(), text.splitlines()))
            num_origins = int(demand_info[0][0])
            od_ptr = list(map(lambda x: [x[0], int(x[1])], demand_info[1:num_origins + 1]))
            des_info = list(map(lambda x: [x[0], float(x[1])*scaling_coeff], demand_info[num_origins + 2:]))

            for i, [origin_i, des_start] in enumerate(od_ptr):
                if i < len(od_ptr) - 1:
                    des_end = od_ptr[i + 1][1] - 1
                else:
                    des_end = len(des_info)

                self.od_demand[origin_i] = dict()
                for des_j in range(des_start - 1, des_end):
                    self.od_demand[origin_i][des_info[des_j][0]] = des_info[des_j][1]
        return self.od_demand

    def expand_bottleneck_capacity(self, scaling_rate, bottleneck_threshold):
        bottleneck_ls = []
        for o in self.graph.adj.keys():
            for d in self.graph[o].keys():
                if self.graph[o][d]['weight'] / self.graph[o][d]['t0'] > bottleneck_threshold:
                    bottleneck_ls.append([o, d])
                    self.graph[o][d]['capacity'] *= scaling_rate
                    self.graph[o][d]['x'] = 0.0
                    self.graph[o][d]['y'] = 0.0
                    self.graph[o][d]['d'] = self.graph[o][d]['y'] - self.graph[o][d]['x']
                    self.graph[o][d]['weight'] = self.graph[o][d]['t0']
        bottleneck_df = pd.DataFrame(bottleneck_ls, columns=['origin', 'end'])
        return bottleneck_df

    def update_link_cost(self, o, d, performance_func='bpr'):
        self.graph[o][d]['weight'] = getattr(self, performance_func)(x=self.graph[o][d]['x'],
                                                                     t0=self.graph[o][d]['t0'],
                                                                     ca=self.graph[o][d]['capacity'])

    def update_x(self, o, d, alpha_hat):
        self.graph[o][d]['x_old'] = self.graph[o][d]['x']
        self.graph[o][d]['x'] += alpha_hat * self.graph[o][d]['d']

    def search_direction(self, o, d):
        self.graph[o][d]['d'] = self.graph[o][d]['y'] - self.graph[o][d]['x']

    def AoN(self, o, d, var_for_update):
        route_nodes, min_cost = label_correcting(self.graph, o, d, Stack())
        for num_node in range(len(route_nodes) - 1):
            self.graph[route_nodes[num_node]][route_nodes[num_node + 1]][var_for_update] += self.od_demand[o][d]

    @staticmethod
    def bpr(x, t0, ca, alpha=0.15, beta=4):
        return t0 * (1 + alpha * pow(x/(0.9 * ca), beta))

    @staticmethod
    def davidson(x, t0, ca, j=0.25):
        return t0 * abs(1+j*(x/(ca-x)))

    def get_objective(self, alpha_hat, flow_pattern, performance_func='bpr'):
        obj_value = 0.0
        for o in self.graph.adj.keys():
            for tail in self.graph[o].keys():
                upper_limit = self.graph[o][tail]['x'] + alpha_hat * self.graph[o][tail]['d']
                if flow_pattern == 'UE':
                    f_x = integrate.quad(getattr(self, performance_func), 0, upper_limit,
                                         args=(self.graph[o][tail]['t0'], self.graph[o][tail]['capacity']))
                    obj_value += f_x[0]

                elif flow_pattern == 'SO':
                    f_x = getattr(self, performance_func)(upper_limit, self.graph[o][tail]['t0'],
                                                          self.graph[o][tail]['capacity']) * upper_limit
                    obj_value += f_x
        return obj_value

    def golden_section(self, a, b, num_iter, flow_pattern, gr=(np.sqrt(5)-1)/2.0, performance_func='bpr'):
        for i in range(num_iter):
            l_alpha = a + (1-gr)*(b-a)
            r_alpha = a + gr * (b-a)

            lhs = self.get_objective(l_alpha, flow_pattern, performance_func=performance_func)
            rhs = self.get_objective(r_alpha, flow_pattern, performance_func=performance_func)

            if lhs < rhs:
                b = r_alpha
            else:
                a = l_alpha
        alpha = (l_alpha + r_alpha)/2.0
        return alpha

    def eval_convergence(self):
        sum_x_old = 0.0
        sum_delta_x = 0.0
        system_cost = 0.0
        for o in self.graph.adj.keys():
            for tail in self.graph[o].keys():
                sum_x_old += self.graph[o][tail]['x_old']
                system_cost += self.graph[o][tail]['x_old'] * self.graph[o][tail]['weight']
                sum_delta_x += np.square(self.graph[o][tail]['x'] - self.graph[o][tail]['x_old'])
        conv = np.sqrt(sum_delta_x)/sum_x_old
        return conv, sum_x_old, system_cost

    def traffic_optimization(self, flow_pattern, conv_threshold, num_iter, performance_func='bpr'):

        run_flag = True
        conv_ls = []

        # initialization of x_i
        for o_i in self.od_demand.keys():
            for d_i in self.od_demand[o_i].keys():
                self.AoN(o_i, d_i, var_for_update='x')

        while run_flag is True:
            # update link-cost
            for o_i in self.graph.adj.keys():
                for tail_i in self.graph[o_i].keys():
                    self.update_link_cost(o_i, tail_i)

            # update y_i
            for o_i in self.od_demand.keys():
                for d_i in self.od_demand[o_i].keys():
                    self.AoN(o_i, d_i, var_for_update='y')

            # search direction d_i
            for o_i in self.graph.adj.keys():
                for tail_i in self.graph[o_i].keys():
                    self.search_direction(o_i, tail_i)

            # line-search
            alpha_optim = self.golden_section(a=0.0, b=1.0, flow_pattern=flow_pattern,
                                              num_iter=num_iter,
                                              performance_func=performance_func)

            # update x
            for o_i in self.graph.adj.keys():
                for tail_i in self.graph[o_i].keys():
                    self.update_x(o_i, tail_i, alpha_hat=alpha_optim)

            conv_res = self.eval_convergence()
            conv_ls.append(conv_res)
            # print('convergence test result is: ', conv_ls)
            if conv_res[0] > conv_threshold:
                for o_i in self.graph.adj.keys():
                    for tail_i in self.graph[o_i].keys():
                        # set y back to zero
                        self.graph[o_i][tail_i]['y'] = 0.0
                continue
            else:
                run_flag = False
                break
        self.conv_df = pd.DataFrame(conv_ls, columns=['conv_value', 'Total Link Flow', 'Total Travel Cost'])


def plot_trend(ue_df, so_df, vars_for_plot):
    x_uplimit = max(ue_df.shape[0], so_df.shape[0])
    for var in vars_for_plot:
        plt.plot(np.arange(1, ue_df.shape[0] + 1), ue_df[var], 'ro-', label='%s of UE' % var)
        plt.plot(np.arange(1, so_df.shape[0] + 1), so_df[var], 'go-', label='%s of SO' % var)
        plt.title('%s vs Iteration Number' % var)
        plt.xlabel('Number of Iteration')
        plt.ylabel('%s' % var)
        plt.xlim(1, x_uplimit)
        plt.xticks(np.arange(1, x_uplimit, 1.0))
        plt.legend()
        plt.show()


def plot_scale_trend(cost_df):

    plt.plot(cost_df['Scaling Coefficient'], cost_df['Total Travel Cost of UE'], 'ro-', label='Total Travel Cost for UE')
    plt.plot(cost_df['Scaling Coefficient'], cost_df['Total Travel Cost of SO'], 'go-', label='Total Travel Cost for SO')
    plt.title('Total Travel Time vs Scaling Coefficient Value')
    plt.xlabel('Coefficient value')
    plt.ylabel('Total Travel Time')
    plt.xlim(min(cost_df['Scaling Coefficient']), max(cost_df['Scaling Coefficient']))
    plt.xticks(cost_df['Scaling Coefficient'])
    plt.legend()
    plt.show()


if __name__ == "__main__":
    # Initialization of traffic conditions
    traff_optim = TrafficOptim()
    traff_optim.load_graph("anaheim.xlsx", ['start', 'end', 'Free Flow Travel Time (min)', 'Capacity'])
    traff_optim.load_demand("fort2.txt")

    conv_threshold = pow(10, -3)
    num_iter = 10  # the line-search iteration number

    # comparison of traffic assignment optimization using UE and SO ===========================
    traff_ue = copy.deepcopy(traff_optim)
    traff_so = copy.deepcopy(traff_optim)

    traff_so.traffic_optimization(flow_pattern='SO', conv_threshold=conv_threshold, num_iter=num_iter)
    traff_ue.traffic_optimization(flow_pattern='UE', conv_threshold=conv_threshold, num_iter=num_iter)

    plot_trend(traff_ue.conv_df, traff_so.conv_df, vars_for_plot=['Total Link Flow', 'Total Travel Cost'])

    # Identify and remove bottlenecks =======================================
    bneck_threshold = 1.15  # the actual link cost is larger than 1.15 times of the free flow travel cost
    scaling_rate = 1.5
    traff_bneck = copy.deepcopy(traff_optim)

    traff_bneck.traffic_optimization(flow_pattern='UE', conv_threshold=conv_threshold, num_iter=num_iter)
    bneck_df = traff_bneck.expand_bottleneck_capacity(scaling_rate=scaling_rate, bottleneck_threshold=bneck_threshold)

    print('The number of identified bottlenecks:', bneck_df.shape[0])
    print('Ten examples of identified bottlenecks:')
    print(bneck_df.head(10))

    # Comparison of scaling parameters =======================================
    scaling_coeff_ls = [0.25, 0.5, 0.75, 1.25, 1.5, 2.0, 3.0, 5.0]
    cost_scale_ls = []
    for coeff in scaling_coeff_ls:
        traff_optim.load_demand("fort2.txt", scaling_coeff=coeff)
        traff_ue_scale = copy.deepcopy(traff_optim)
        traff_so_scale = copy.deepcopy(traff_optim)

        traff_ue_scale.traffic_optimization(flow_pattern='UE', conv_threshold=conv_threshold, num_iter=num_iter)
        traff_so_scale.traffic_optimization(flow_pattern='SO', conv_threshold=conv_threshold, num_iter=num_iter)

        cost_scale_ls.append([coeff,
                             traff_ue_scale.conv_df.iloc[-1]['Total Travel Cost'],
                             traff_so_scale.conv_df.iloc[-1]['Total Travel Cost']])
    cost_scale_df = pd.DataFrame(cost_scale_ls, columns=['Scaling Coefficient', 'Total Travel Cost of UE',
                                                         'Total Travel Cost of SO'])
    plot_scale_trend(cost_scale_df)
