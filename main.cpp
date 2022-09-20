#include "graph.h"
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <cmath>

using namespace std;

clock_t start_time;
graph G;
int n_original;
ofstream f_log;
bool parallel_compression;
int n_tries;
double p_tries;
bool check3;
vector <double> tmp_w;
void get_start_time()
{
    start_time = clock();
}

double get_running_time()
{
    return double(clock() - start_time) / (CLOCKS_PER_SEC);
}

void print_graph_test()
{
    int cnt;
    for (int i = 0; i < G.n_nodes; ++i)
    {
        // f_log << i << ": ";
        // for (int j: G.nodes[i].ID)
        //     f_log << j << " ";
        // f_log << endl;
        cnt = 0;
        f_log << i << " " << G.nodes[i].s<< endl;
        for (
            
            auto j: G.nodes[i].edges)
        {
            // f_log << i << " " << j.neighbor << " " << j.w << " " << j.pos << " " << cnt << endl;
            f_log << i << " " << j.neighbor << " " << j.w << endl;            
            ++cnt;
        }
        f_log << "--------------------------------------\n";
    }
    // for (int i = 0; i < n_original; ++i)
    //     f_log << G.is_flipped[i] << " ";
    // f_log << endl;
    f_log << "--------------------------------------\n";
}

void check_graph()
{
    for (int i = 0; i < G.n_nodes; ++i)
    {        
        for (auto j: G.nodes[i].edges)
        {
            if (j.neighbor >= G.n_nodes) 
            {
                cerr << "WRONG!!!!\n";            
                return;
            }
        }        
    }
}

bool fast_find_NS()
{
    bool compressed = false;
    vector <int> g;
    g.resize(2);
    int u,v,z;
    double w,w1,w2;
    vector <vector <int> > links;
    links.resize(G.n_nodes);
    vector <int> color(G.n_nodes,0);
    vector <int> ng;    
    vector <int> weak_ng;
    bool is_flipped_weak;
    bool antipolar;
    double _nsi;
    // cerr << G.n_nodes << endl;
    f_log << "Find NGs with k = 2\n";
    clock_t time_tmp;
    time_tmp = clock();
    for (u = 0; u < G.n_nodes; ++u)
        if (G.nodes[u].edges.size() == 0)
        {
            G.swap_nodes(u,G.n_nodes-1);
            G.delete_node();
        }
    for (u = 0; u < G.n_nodes; ++u)
    {      
        if (G.nodes[u].l < G.l) continue;  
        for (auto e: G.nodes[u].edges)
        {
            v = e.neighbor;
            if (u >= v && G.nodes[v].l == G.l) continue;
            w = e.w;
            _nsi = 2*abs(w) - min(G.nodes[u].s,G.nodes[v].s);
            if (_nsi > 0)
            {
                if (w < 0)
                {
                    links[u].push_back(G.n_nodes + v);
                    links[v].push_back(G.n_nodes + u);
                    // cout << "antipolar :" << u << " " << v << endl;
                }
                else
                {
                    links[u].push_back(v);
                    links[v].push_back(u);
                    // cout << "Non-separate :" << u << " " << v << endl;
                }
            }
            if ((_nsi == 0) && (weak_ng.size() == 0))
            {
                G.nodes[u].l = G.l;
                weak_ng.push_back(v);
                weak_ng.push_back(u);  
                if (w < 0)
                    is_flipped_weak = true;
                else is_flipped_weak = false;
            } 
        }
    }
    ++G.l;
    f_log << "Finding NGs in " << double(clock() - time_tmp) / CLOCKS_PER_SEC << "s\n";
    time_tmp = clock();
    // cerr << G.n_nodes << endl;
    for (u = 0; u < G.n_nodes; ++u)
    {
        G.pos[u] = u;
        G.rpos[u] = u;
    }

    int n_tmp = G.n_nodes;
    for (u = 0; u < G.n_nodes; ++u)
    {
        if ((color[u] == 0) && (links[u].size() > 0))
        {
            color[u] = 1;
            ng.push_back(u);
            int x,y;
            int cnt = 0;
            // cout << "Flip: ";
            while (cnt < ng.size())
            {
                x = ng[cnt];
                for(int v: links[x])
                {
                    if (color[v % n_tmp] != 0) continue;
                    if (v >= n_tmp)
                    {
                        v -= n_tmp;                        
                        color[v] = -color[x];
                    }
                    else color[v] = color[x];
                    ng.push_back(v);
                    if (color[v] == -1)
                    {
                        G.flip(G.pos[v]);
                        // cout << G.pos[v] << " ";
                    }
                }
                ++cnt;
            }
            // cout << endl;
            // cout << "Pos of 56: " << G.pos[56] << endl;
            // cout << "NG (" << G.n_nodes << "): ";
            // for (int i = 0; i < int(ng.size()); ++i)
            //     cout << ng[i] << " ";            
            // cout << endl;         
            for (int i = 0; i < int(ng.size()); ++i)
                ng[i] = G.pos[ng[i]];
            // cout << "NG (" << G.n_nodes << "): ";
            // for (int i = 0; i < int(ng.size()); ++i)
            //     cout << ng[i] << " ";            
            // cout << endl;         
            G.merge_multiple_nodes(ng);
            // check_graph();
            // print_graph_test();
            ng.clear();
            compressed = true;
        }
    }
    // cerr << "DONE strong" << endl;
    // return false;
    if (compressed == true)
    {
        f_log << "Merging NGs in " << double(clock() - time_tmp) / CLOCKS_PER_SEC << "s\n";
        return true;
    }
    if (weak_ng.size() > 0)
    {
        // cerr << "Start weak" << endl;
        // cerr << weak_ng.size() << endl; 
        if (is_flipped_weak) G.flip(weak_ng[1]);
        G.merge_multiple_nodes(weak_ng);
        // cerr << "DONE weak" << endl;
        return true;
    }
    return false;
}

pair < double,int> _nsi3(double w1, double w2, double w3, double a1, double a2, double a3, int u, int v, int z)
{
    double _w1 = w1, _w2 = w2, _w3 = w3;
    double n_w1, n_w2, n_w3;
    int c = 0;
    if (w1 < 0)
    {
        ++c;
        _w1 = -w1;
    }
    if (w2 < 0)
    {
        ++c;
        _w2 = -w2;
    }
    if (w3 < 0)
    {
        ++c;
        _w3 = -w3;
    }
    a1 = a1 - _w1 - _w2;
    a2 = a2 - _w1 - _w3;
    a3 = a3 - _w3 - _w2;
    n_w1 = _w1;
    n_w2 = _w2;
    n_w3 = _w3;
    if (c % 2 == 1)
    {
        if (_w1 < _w2)
        {
            if (_w1 < _w3)
                n_w1 = -_w1;
            else n_w3 = -_w3;
        }
        else
        {
            if (_w2 < _w3)
                n_w2 = -_w2;
            else n_w3 = -_w3;
        }
    }
    pair<double,int> res;
    
    res.first = min(n_w1+n_w2-min(a1,a2+a3),n_w1+n_w3-min(a2,a1+a3));
    res.first = min(res.first, n_w2+n_w3-min(a3,a1+a2));
    res.second = -1;
    if (res.first >=0)
    {
        if (w1 ==0)
        {
            if (w2 != n_w2 && w3 != n_w3) res.second = z;
            else if (w2 != n_w2) res.second = u;
            else if (w3 != n_w3) res.second = v;
        } else
        if (w2 ==0)
        {
            if (w1 != n_w1 && w3 != n_w3) res.second = v;
            else if (w1 != n_w1) res.second = u;
            else if (w3 != n_w3) res.second = z;
        } else
        if (w3 ==0)
        {
            if (w2 != n_w2 && w1 != n_w1) res.second = u;
            else if (w1 != n_w1) res.second = v;
            else if (w2 != n_w2) res.second = z;
        } else
        if (n_w1!=w1) 
        {
            if (n_w2!=w2) res.second = u;
            else res.second = v;
        }
        else
        {
            if (n_w2!=w2) res.second = z;
        }
        // f_log << n_w1+n_w2-min(a1,a2+a3) << " " << n_w1+n_w3-min(a2,a1+a3) << " " << n_w2+n_w3-min(a3,a1+a2) << endl;
        // f_log << res.first << " " << res.second << endl;
        // f_log << u << " " << v << " " << z << endl;
        // f_log << w1 << " " << w2 << " " << w3 << endl;
        // f_log << _w1 << " " << _w2 << " " << _w3 << endl;
        // f_log << n_w1 << " " << n_w2 << " " << n_w3 << endl;
        // f_log << "______________________" << endl;
    }    
    return res;
}

void find3(int u, int v, double w1, vector <vector <int> > &links)
{
    int z;
    pair<double,int> res;;
    // double w1,w2;
    for (auto e: G.nodes[u].edges)
    {
        z = e.neighbor;
        tmp_w[z] = e.w;
    }    
    for (auto e: G.nodes[v].edges)
    {
        z = e.neighbor;
        if (z == u) continue;
        res = _nsi3(w1,tmp_w[z],e.w,G.nodes[u].s, G.nodes[v].s, G.nodes[z].s,u,v,z);
        if (res.first > 0)
        {
            if (res.second == -1)
            {
                links[u].push_back(v);
                links[v].push_back(u);
                links[u].push_back(z);
                links[z].push_back(u);
            } 
            else if (res.second == u)
            {
                links[u].push_back(G.n_nodes + v);
                links[v].push_back(G.n_nodes + u);
                links[u].push_back(G.n_nodes + z);
                links[z].push_back(G.n_nodes + u);
            }
            else if (res.second == v)
            {
                links[u].push_back(G.n_nodes + v);
                links[v].push_back(G.n_nodes + u);
                links[u].push_back(z);
                links[z].push_back(u);
            }
            else if (res.second == z)
            {
                links[u].push_back(v);
                links[v].push_back(u);
                links[u].push_back(G.n_nodes + z);
                links[z].push_back(G.n_nodes + u);
            }
        }
    }
    for (auto e: G.nodes[u].edges)
    {
        z = e.neighbor;
        tmp_w[z] = 0;
    }
}

bool slow_find_NS()
{
    // two nodes
    bool compressed = false;
    vector <int> g;
    g. resize(2);
    int u,v,z;
    double w,w1,w2;
    vector <vector <int> > links;
    links.resize(G.n_nodes);
    vector <int> color(G.n_nodes,0);
    vector <int> ng;    
    vector <int> weak_ng;
    bool is_flipped_weak;
    bool antipolar;
    double _nsi;
    // cerr << G.n_nodes << endl;
    f_log << "Find NGs with k = 2\n";
    clock_t time_tmp;
    time_tmp = clock();
    vector <pair<double,int> > sc;
    int cnt = 0;
    vector <int> uu, vv;
    vector <double> ww;
    n_tries = 0;     
    for (u = 0; u < G.n_nodes; ++u)
        if (!G.check_slow || G.nodes[u].l >= G.l-1)
        { 
            ++n_tries;       
            for (auto e: G.nodes[u].edges)
            {
                v = e.neighbor;
                if (u >= v) continue;
                w = e.w;
                _nsi = 2*abs(w) - min(G.nodes[u].s,G.nodes[v].s);
                sc.push_back(make_pair(_nsi,cnt));
                ++cnt;
                uu.push_back(u);
                vv.push_back(v);
                ww.push_back(w);
            }
        }    
    G.check_slow = true;
    sort(sc.rbegin(),sc.rend());
    n_tries = min(int(n_tries*p_tries),int(sc.size()));
    // cout << n_tries << " " << p_tries << endl;
    for (int i = 0; i < n_tries;++i)
    {
        u = uu[sc[i].second];
        v = vv[sc[i].second];
        w = ww[sc[i].second];
        if (check3)
        {
            find3(u,v,w,links);
            find3(v,u,w,links);
        }
        antipolar = false;
        if (w < 0)
        {
            G.flip(u);
            antipolar = true;
        }
        g[0] = u;
        g[1] = v;
        _nsi = G.NSI(g);
        if (_nsi > 0)
        {
            if (antipolar)
            {
                links[u].push_back(G.n_nodes + v);
                links[v].push_back(G.n_nodes + u);
                // cout << "antipolar :" << u << " " << v << endl;
            }
            else
            {
                links[u].push_back(v);
                links[v].push_back(u);
                // cout << "Non-separate :" << u << " " << v << endl;
            }
        }
        if ((_nsi == 0) && (weak_ng.size() == 0))
        {
            weak_ng.push_back(v);
            weak_ng.push_back(u);  
            is_flipped_weak = antipolar;              
        }
        if (antipolar)
        {
            G.flip(u);
        }

    }
    f_log << "Finding NGs in " << double(clock() - time_tmp) / CLOCKS_PER_SEC << "s\n";
    time_tmp = clock();
    // cerr << G.n_nodes << endl;
    for (u = 0; u < G.n_nodes; ++u)
    {
        G.pos[u] = u;
        G.rpos[u] = u;
    }
    // for (u = 0; u < G.n_nodes; ++u)
    // {
    //     f_log << "links\n";
    //     f_log << "------------------\n";
    //     f_log << u << endl;
    //     for (int ii: links[u])
    //         f_log << v << " ";
    //     f_log << "\n------------------\n";
    // }
    int n_tmp = G.n_nodes;
    for (u = 0; u < G.n_nodes; ++u)
    {
        if ((color[u] == 0) && (links[u].size() > 0))
        {
            color[u] = 1;
            ng.push_back(u);
            int x,y;
            int cnt = 0;
            // cout << "Flip: ";
            while (cnt < ng.size())
            {
                x = ng[cnt];
                for(int v: links[x])
                {
                    if (color[v % n_tmp] != 0) continue;
                    if (v >= n_tmp)
                    {
                        v -= n_tmp;                        
                        color[v] = -color[x];
                    }
                    else color[v] = color[x];
                    ng.push_back(v);
                    if (color[v] == -1)
                    {
                        G.flip(G.pos[v]);
                        // cout << G.pos[v] << " ";
                    }
                }
                ++cnt;
            }
            // cout << endl;
            // cout << "Pos of 56: " << G.pos[56] << endl;
            // cout << "NG (" << G.n_nodes << "): ";
            // for (int i = 0; i < int(ng.size()); ++i)
            //     cout << ng[i] << " ";            
            // cout << endl;         
            for (int i = 0; i < int(ng.size()); ++i)
                ng[i] = G.pos[ng[i]];
            // cout << "NG (" << G.n_nodes << "): ";
            // for (int i = 0; i < int(ng.size()); ++i)
            //     cout << ng[i] << " ";            
            // cout << endl;         
            G.merge_multiple_nodes(ng);
            // check_graph();
            // print_graph_test();
            ng.clear();
            compressed = true;
        }
    }
    // cerr << "DONE strong" << endl;
    // return false;
    if (compressed == true)
    {
        f_log << "Merging NGs in " << double(clock() - time_tmp) / CLOCKS_PER_SEC << "s\n";
        return true;
    }
    if (weak_ng.size() > 0)
    {
        // cerr << "Start weak" << endl;
        // cerr << weak_ng.size() << endl; 
        if (is_flipped_weak) G.flip(weak_ng[1]);
        G.merge_multiple_nodes(weak_ng);
        // cerr << "DONE weak" << endl;
        return true;
    }


    return false;
}

bool find_NS(int k)
{
    // two nodes
    bool compressed = false;
    vector <int> g;
    g. resize(2);
    int u,v,z;
    double w,w1,w2;
    vector <vector <int> > links;
    links.resize(G.n_nodes);
    vector <int> color(G.n_nodes,0);
    vector <int> ng;    
    vector <int> weak_ng;
    bool is_flipped_weak;
    bool antipolar;
    double _nsi;
    // cerr << G.n_nodes << endl;
    f_log << "Find NGs with k = 2\n";
    clock_t time_tmp;
    time_tmp = clock();
    for (u = 0; u < G.n_nodes; ++u)
    {        
        for (auto e: G.nodes[u].edges)
        {
            v = e.neighbor;
            if (u >= v) continue;
            w = e.w;
            // cerr << u << " " << v << " " << G.n_nodes << endl;
            antipolar = false;
            if (w < 0)
            {
                G.flip(u);
                antipolar = true;
            }
            g[0] = u;
            g[1] = v;
            _nsi = G.NSI(g);
            if (!parallel_compression)
            {
                if (_nsi >= 0)
                {
                    f_log << "Finding NGs in " << double(clock() - time_tmp) / CLOCKS_PER_SEC << "s\n";
                    time_tmp = clock();
                    G.merge_multiple_nodes(g);
                    f_log << "Merging NGs in " << double(clock() - time_tmp) / CLOCKS_PER_SEC << "s\n";
                    return true;
                }
            }
            if (_nsi > 0)
            {
                if (antipolar)
                {
                    links[u].push_back(G.n_nodes + v);
                    links[v].push_back(G.n_nodes + u);
                    // cout << "antipolar :" << u << " " << v << endl;
                }
                else
                {
                    links[u].push_back(v);
                    links[v].push_back(u);
                    // cout << "Non-separate :" << u << " " << v << endl;
                }
            }
            if ((_nsi == 0) && (weak_ng.size() == 0))
            {
                weak_ng.push_back(v);
                weak_ng.push_back(u);  
                is_flipped_weak = antipolar;              
            }
            if (antipolar)
            {
                G.flip(u);
            }
        }
    }
    f_log << "Finding NGs in " << double(clock() - time_tmp) / CLOCKS_PER_SEC << "s\n";
    time_tmp = clock();
    // cerr << G.n_nodes << endl;
    for (u = 0; u < G.n_nodes; ++u)
    {
        G.pos[u] = u;
        G.rpos[u] = u;
    }

    int n_tmp = G.n_nodes;
    for (u = 0; u < G.n_nodes; ++u)
    {
        if ((color[u] == 0) && (links[u].size() > 0))
        {
            color[u] = 1;
            ng.push_back(u);
            int x,y;
            int cnt = 0;
            // cout << "Flip: ";
            while (cnt < ng.size())
            {
                x = ng[cnt];
                for(int v: links[x])
                {
                    if (color[v % n_tmp] != 0) continue;
                    if (v >= n_tmp)
                    {
                        v -= n_tmp;                        
                        color[v] = -color[x];
                    }
                    else color[v] = color[x];
                    ng.push_back(v);
                    if (color[v] == -1)
                    {
                        G.flip(G.pos[v]);
                        // cout << G.pos[v] << " ";
                    }
                }
                ++cnt;
            }
            // cout << endl;
            // cout << "Pos of 56: " << G.pos[56] << endl;
            // cout << "NG (" << G.n_nodes << "): ";
            // for (int i = 0; i < int(ng.size()); ++i)
            //     cout << ng[i] << " ";            
            // cout << endl;         
            for (int i = 0; i < int(ng.size()); ++i)
                ng[i] = G.pos[ng[i]];
            // cout << "NG (" << G.n_nodes << "): ";
            // for (int i = 0; i < int(ng.size()); ++i)
            //     cout << ng[i] << " ";            
            // cout << endl;         
            G.merge_multiple_nodes(ng);
            // check_graph();
            // print_graph_test();
            ng.clear();
            compressed = true;
        }
    }
    // cerr << "DONE strong" << endl;
    // return false;
    if (compressed == true)
    {
        f_log << "Merging NGs in " << double(clock() - time_tmp) / CLOCKS_PER_SEC << "s\n";
        return true;
    }
    if (weak_ng.size() > 0)
    {
        // cerr << "Start weak" << endl;
        // cerr << weak_ng.size() << endl; 
        if (is_flipped_weak) G.flip(weak_ng[1]);
        G.merge_multiple_nodes(weak_ng);
        // cerr << "DONE weak" << endl;
        return true;
    }
    // three nodes
    if (k >= 3)
    {
        g.resize(3);    
        for (u = 0; u < G.n_nodes; ++u)
        {
            int n = G.nodes[u].edges.size();
            for (int i = 0; i < n; ++i)
            {
                v = G.nodes[u].edges[i].neighbor;
                if (u >= v) continue;
                for (int j = i+1; j < n; ++j)
                {
                    z = G.nodes[u].edges[j].neighbor;
                    if (u >= z) continue;
                    g[0] = u;
                    g[1] = v;
                    g[2] = z;
                    _nsi = G.NSI(g);
                    if (_nsi >= 0)
                    {
                        G.merge_multiple_nodes(g);
                        return true;
                    }
                }                
            }            
        }
    }    

    return false;
}

bool find_probe()
{
    bool compressed = false;
    int u , x, y;
    vector <int> adj;
    vector <int> g(2);
    vector <int> ng;
    vector <int> ngs(G.n_nodes,-1);
    node a,b;
    bool antipolar;
    double _nsi;
    vector <int> pp;
    bool check_flip;
    for (u = 0; u < G.n_nodes; ++u)
    {
        adj.clear();
        // cout << "check :" << " " << u << " " << G.n_nodes<< endl;
        for (auto e: G.nodes[u].edges)
        {            
            if (u < e.neighbor) adj.push_back(e.neighbor);
        }
        for (int v : adj)
        {
            // cout << "check :" << " " << u << " " << v << " "<< G.n_nodes<< endl;
            G.merge_nodes(u,v,a,b);
            x = G.n_nodes - 1;
            for (auto e: G.nodes[x].edges)
            {
                y = e.neighbor;
                antipolar = false;
                if (e.w < 0)
                {
                    G.flip(x);
                    antipolar = true;
                }
                g[0] = x;
                g[1] = y;
                _nsi = G.NSI(g);
                if (antipolar == true)
                    G.flip(x);
                if (_nsi >= 0)
                {
                    if (antipolar) ngs[y] = 1;
                    else ngs[y] = 0;
                    pp.push_back(y);
                }
            }
            G.restore(a,b);

            G.flip(G.n_nodes-2);
            G.merge_nodes(G.n_nodes-1,G.n_nodes-2,a,b);
            for (auto e: G.nodes[x].edges)
            {
                y = e.neighbor;
                if (ngs[y] == -1) continue;
                // cout << "check :" << u <<  " " << v << " " << x << " " << y << endl;
                antipolar = false;
                if (e.w < 0)
                {
                    G.flip(x);
                    antipolar = true;
                }
                g[0] = x;
                g[1] = y;
                _nsi = G.NSI(g);
                if (antipolar == true)
                    G.flip(x);
                if (_nsi >= 0)
                {
                    ng.resize(2);
                    if (ngs[y] == 1)
                    {
                        check_flip = true;
                        if (antipolar) 
                        {
                            ng[0] = G.n_nodes;
                            ng[1] = y;
                        }
                        else
                        {
                            ng[0] = G.n_nodes-1;
                            ng[1] = y;   
                        }
                    }
                    else
                    {
                        check_flip = false;
                        if (antipolar) 
                        {
                            ng[0] = G.n_nodes-1;
                            ng[1] = y;
                        }
                        else
                        {
                            ng[0] = G.n_nodes;
                            ng[1] = y;   
                        }
                    }
                }
                // if (ng.size() > 0) break;
            }
            G.restore(a,b);
            G.flip(G.n_nodes-2);
            for (int z: pp)
                ngs[z] = -1;
            pp.clear();
            if (ng.size() > 0) break;
            G.restore_order(u,v);
        }
        if (ng.size() > 0) break;
    }
    if (ng.size() > 0)
    {
        // for (int i = 0; i < int(ng.size()); ++i)
        // {
        //     if (ng[i] < 0)
        //     {
        //         ng[i] = abs(ng[i]);
        //         G.flip(ng[i]);
        //     }
        // }
        // cout <<  ng[0] << " " << ng[1] << " " << check_flip << endl;
        if (check_flip) G.flip(ng[1]);
        G.merge_multiple_nodes(ng);
        // print_graph_test();
        return true;
    }
    return false;
}


void compress(int k, bool probe)
{
    bool compressed = true;
    int cnt = 0;
    clock_t time_tmp;
    int tmp;
    while (compressed && (G.n_nodes > 1) )
    {
        ++cnt;
        f_log << "Iteration " << cnt << " start with " << G.n_nodes << " nodes\n";
        print_graph_test();
        tmp = G.n_nodes;
        time_tmp = clock();
        compressed = fast_find_NS();
        // compressed = find_NS(k);
        print_graph_test();
        if (compressed == false || G.check_slow)
        {   
            compressed = slow_find_NS();             
        }
        f_log << "Iteration " << cnt << " compress " << tmp - G.n_nodes << " nodes in " << double(clock() - time_tmp) / CLOCKS_PER_SEC << "s\n";
    }
}

void print_graph(string out_file)
{
    ofstream f1,f2;
    string s1 = out_file+"_compress";
    string s2 = out_file+"_flip";
    f1.open(s1);
    f2.open(s2);
    int m = 0;
    f1 << G.n_nodes << " ";    
    for (int i = 0; i < G.n_nodes; ++i)
        for (auto j : G.nodes[i].edges)
            if ((i < j.neighbor) && (fabs(j.w)>1e-8))
                m = m + 1;
    f1 << m << endl;
    for (int i = 0; i < G.n_nodes; ++i)
        for (auto j : G.nodes[i].edges)
            if ((i < j.neighbor) && (fabs(j.w)>1e-8))
                f1 << i << " " << j.neighbor << " " << j.w << endl;
    f2 << G.n_nodes << " " << n_original << endl;
    for (int i = 0; i < G.n_nodes; ++i)
    {
        f2 << G.nodes[i].ID.size() << " ";
        for (int j: G.nodes[i].ID)
            f2 << j << " ";
        f2 << endl;
    }   
    for (int i = 0; i < n_original; ++i)
        f2 << int(G.is_flipped[i]) << " ";
    f2 << endl;
    f1.close();
    f2.close();    
}


int main(int argc, char ** argv)
{
    // pair<double,int> res;
    // double w1,w2,w3,a1,a2,a3;
    // cin >> w1 >> w2 >> w3 >> a1 >> a2 >> a3;
    // res = _nsi3(w1,w2,w3,a1,a2,a3,1,2,3);
    // return 0;    
    if (argc < 2) {
        cout << argv[0]<<" <input-file in .net format> <output-file> alpha"<<endl;
        return 0;
    }
    cout << argv[1] <<","; 
    string log_file = argv[2];
    log_file = log_file + "_log";
    f_log.open(log_file);
    G.load_graph(argv[1]);
    tmp_w.resize(G.n_nodes,0);
    n_original = G.n_nodes;
    bool probe = false;
    p_tries = 0.25;
    if (argc > 3) 
        p_tries = atof(argv[3]);
    n_tries = int(p_tries*G.n_nodes);
    check3 = true;
    //if (argc >4 && atoi(argv[4]) > 0) check3 = true;
    // int x = atoi(argv[3]);
    // if (x > 0) probe = true;
    int k = 2;
    // if (argc > 4)
    // {
    //     k = atoi(argv[4]);        
    // }
    parallel_compression = true;
    // if (argc > 5)
    // {
    //     if (atoi(argv[5]) == 0)
    //         parallel_compression = false;
    // }
    node a,b;
    // a = G.nodes[3];
    // G.delete_node();
    // b = G.nodes[2];  
    // G.delete_node();
    // print_graph_test();
    // G.add_node(b);
    // G.add_node(a);
    // print_graph_test();
    // G.swap_nodes(0,3);
    // print_graph_test();
    // G.merge_nodes();
    // print_graph_test();
    // G.merge_nodes(0,1,a,b);    
    // G.restore(a,b);
    // G.flip(2);
    // G.merge_nodes(2,3,a,b);
    // G.restore(a,b);
    // G.flip(2);
    // G.restore_order(0,1);
    // print_graph_test();
    // return 0;
    // G.restore(0,1,a,b);
    // print_graph_test();
    // return 0;
    // G.flip(0);
    // print_graph_test();
    // G.flip(1);
    // print_graph_test();
    // vector<int> g;
    // int n;
    // cin >> n;
    // g.resize(n);
    // for (int i = 0; i < n; ++i)
    //     cin >> g[i]; 
    // cout << G.NSI(g) << endl;
    // G.merge_multiple_nodes(g);
    // print_graph_test();    
    get_start_time();
    compress(k, probe);
    cout << G.n_nodes <<",";
    cout <<get_running_time() << endl;  
    string s = argv[2];
    // s = s + "_" + to_string(probe) + "_" + to_string(k);
    print_graph(s);  
    f_log.close();
    return 0;
}