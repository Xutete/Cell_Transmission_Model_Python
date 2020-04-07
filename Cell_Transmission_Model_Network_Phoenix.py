# -*- coding: utf-8 -*-
"""
Created on Mon Sep  9 17:33:32 2019

@author: Lyy
"""

import pandas as pd
import numpy as np
import re

class Cell(object):
    idcase = {}
    def __init__(self, cellid, linkid, zoneid, time_interval=6, k=0, qmax=2160, kjam=220, vf=60, w=12, length=0.1, updated=False, arr_rate=0, dis_rate=0):
        self.kjam = kjam
        self.cellid = cellid # local address
        self.linkid = linkid # link layer address
        self.zoneid = zoneid # zone layer address
        self.vf = vf # Time interval = length / vf
        self.w = w
        self.cfrom = []
        self.cto = []
        self.k = k # density at time interval t
        self.oldk = k # density at time interval t-1
        self.qmax = qmax
        self.length = length
        self.updated = updated
        self.arr_rate = arr_rate
        self.dis_rate = dis_rate
        self.time_sec = time_interval
        self.time_hour = time_interval / 3600
        self.inflow = 0
        self.outflow = 0
        self.pk = 0.75
        self.pck = 0.25
        if Cell.idcase.get(self.getCompleteAddress()) == None:
            Cell.idcase.setdefault(self.getCompleteAddress(), self)
        else:
            raise Exception("This id has been used by other cell")
        
    def addConnection(self, sink):
        if len(sink.cfrom) == 2 or len(self.cto) == 2:
            raise Exception("Cannot add more connection to cell %s and cell %s" % (self.getCompleteAddress(), sink.getCompleteAddress())) 
            
        if (len(self.cto) and len(sink.cfrom)) and (len(sink.cto) == 2 or len(self.cfrom) == 2):
            raise Exception("Invaild cell connection! A cell cannot connect to merge and diverge cell simultaneously")
            
        self.cto.append(sink) # An instance of cell class is stored, in order to use cto and cfrom as pointer.
        sink.cfrom.append(self)
        
        if len(self.cto) >= 2:
            interSection(self)
        elif len(sink.cfrom) >= 2:
            interSection(sink)
        
    def deleteConnection(self, sink):
        if sink not in self.cto:
            raise Exception("Cell %s is not connected with cell %s" % (self.getCompleteAddress(), sink.getCompleteAddress()))
            
        self.cto.remove(sink)
        sink.cfrom.remove(self)
        
    def getCell(cid):
        return Cell.idcase[cid]
    
    def getFirstCell(linkid):
        newDict = {}
        for key in Cell.idcase:
            if Cell.idcase[key].linkid == linkid:
                newDict[key] = Cell.idcase[key]
                
        return newDict[min(newDict.keys())]
    
    def getLastCell(linkid):
        newDict = {}
        for key in Cell.idcase:
            if Cell.idcase[key].linkid == linkid:
                newDict[key] = Cell.idcase[key]
                
        return newDict[max(newDict.keys())]
    
    def deleteCell(cid):
        poped = Cell.idcase.pop(cid)
        for elem in poped.cto:
            poped.deleteConnection(elem)
        del poped
        
    def getCompleteAddress(self):
        return "%s.%s.%s" % (self.zoneid, self.linkid, self.cellid)
       
    def updateDensity(self): # This method can only be used by normal cell instance.
        if not self.updated:
            self.oldk = self.k
        if len(self.cfrom) == 2: # Merge at here, we need to update density among this cell and two other upstream cells.
            pk = self.pk # probability from upstream normal cell
            pck = 1 - self.pk # probability from upstream merge cell
            for elem in self.cfrom:
                rek = np.min([self.qmax, self.w * (self.kjam - self.oldk)]) * self.time_hour / self.length
                if elem.linkid == self.linkid:
                    sbk = np.min([elem.qmax, elem.vf * elem.oldk]) * elem.time_hour / elem.length
                    prov = elem
                    
                else:
                    sck = np.min([elem.qmax, elem.vf * elem.oldk]) * elem.time_hour / elem.length
                    if not elem.updated:
                        elem.oldk = elem.k
                        
                    merge = elem
            
            try: # In order to cope with situation that provious cell is the first cell (cfrom is empty)
                prov.inflow = np.min([prov.qmax, prov.vf * prov.cfrom[0].oldk, prov.w * (prov.kjam - prov.oldk)]) * prov.time_hour / prov.length
                prov.outflow = np.min([np.median([pk * rek, sbk, rek - sck]), prov.vf * prov.oldk * prov.time_hour / prov.length])
                
            except:
                prov.inflow = np.min([prov.qmax, prov.arr_rate, prov.w * (prov.kjam - prov.oldk)]) * prov.time_hour / prov.length
                prov.outflow = np.min([np.median([pk * rek, sbk, rek - sck]), prov.vf * prov.oldk * prov.time_hour / prov.length])
            
            if len(merge.cfrom):                
                merge.inflow = np.min([merge.qmax, merge.vf * merge.cfrom[0].oldk, merge.w * (merge.kjam - merge.oldk)]) * merge.time_hour / merge.length
                merge.outflow = np.min([np.median([pck * rek, sck, rek - sbk]), merge.vf * merge.oldk * merge.time_hour / merge.length])
            else:
                merge.inflow = np.min([merge.qmax, merge.arr_rate, merge.w * (merge.kjam - merge.oldk)]) * merge.time_hour / merge.length
                merge.outflow = np.min([np.median([pck * rek, sck, rek - sbk]), merge.vf * merge.oldk * merge.time_hour / merge.length])
            
            if len(self.cto):                
                self.inflow = np.min([self.qmax * self.time_hour / self.length, sbk+sck, self.w * (self.kjam - self.oldk) * self.time_hour / self.length])
                self.outflow = np.min([self.cto[0].qmax, self.oldk * self.vf, self.cto[0].w * (self.cto[0].kjam - self.cto[0].oldk)]) * self.time_hour / self.length
            else:
                self.inflow = np.min([self.qmax * self.time_hour / self.length, sbk+sck, self.w * (self.kjam - self.oldk) * self.time_hour / self.length])
                self.outflow = np.min([self.qmax, self.oldk * self.vf, self.dis_rate]) * self.time_hour / self.length
            
            prov.k = prov.oldk + np.max([0, prov.inflow]) - np.max([0, prov.outflow])
            merge.k = merge.oldk + np.max([0, merge.inflow]) - np.max([0, merge.outflow])
            self.k = self.oldk + np.max([0, self.inflow]) - np.max([0, self.outflow])
            
            prov.updated, self.updated, merge.updated = True, True, True
            
        elif len(self.cto) == 2: # Diverge at here
            ptnc = self.pk # Propotion towards to next normal cell
            ptdc = 1 - self.pk # Propotion towards to diverge cell
            for elem in self.cto:
                if elem.linkid == self.linkid:
                    elem.oldk = elem.k
                    next_c = elem
                
                else:
                    if not elem.updated:
                        elem.oldk = elem.k
                        
                    diverge = elem
            
            rck = np.min([next_c.qmax, next_c.w * (next_c.kjam - next_c.oldk)]) * next_c.time_hour / next_c.length # Receive ability of next normal cell
            rek = np.min([diverge.qmax, diverge.w * (diverge.kjam - diverge.oldk)]) * diverge.time_hour / diverge.length
            sbk = np.min([self.qmax, self.vf * self.oldk]) * self.time_hour / self.length
            
            try:# In order to cope with situation that next cell is the last cell (cto is empty)
                next_c.inflow = ptnc * np.min([sbk, rek/ptdc, rck/ptnc])
                next_c.outflow = np.min([next_c.cto[0].qmax, next_c.vf * next_c.oldk, next_c.cto[0].w * (next_c.cto[0].kjam - next_c.cto[0].oldk)]) * next_c.time_hour / next_c.length
            except:
                next_c.inflow = ptnc * np.min([sbk, rek/ptdc, rck/ptnc])
                next_c.outflow = np.min([next_c.qmax, next_c.vf * next_c.oldk, next_c.dis_rate]) * next_c.time_hour / next_c.length
            
            if len(diverge.cto):
                diverge.inflow = ptdc * np.min([sbk, rek/ptdc, rck/ptnc])
                diverge.outflow = np.min([diverge.cto[0].qmax, diverge.oldk * diverge.vf, diverge.cto[0].w * (diverge.cto[0].kjam - diverge.cto[0].oldk)]) * diverge.time_hour / diverge.length
            else:
                diverge.inflow = ptdc * np.min([sbk, rek/ptdc, rck/ptnc])
                diverge.outflow = np.min([diverge.qmax, diverge.oldk * diverge.vf, diverge.dis_rate]) * diverge.time_hour / diverge.length
            
            if len(self.cfrom):
                self.inflow = np.min([self.qmax, self.cfrom[0].oldk * self.vf, self.w * (self.kjam - self.oldk)]) * self.time_hour / self.length
                self.outflow = np.min([sbk, rek/ptdc, rck/ptnc])
            else:
                self.inflow = np.min([self.qmax, self.arr_rate, self.w * (self.kjam - self.oldk)]) * self.time_hour / self.length
                self.outflow = np.min([sbk, rek/ptdc, rck/ptnc])
            
            next_c.k = next_c.oldk + np.max([0, next_c.inflow]) - np.max([0, next_c.outflow])
            diverge.k = diverge.oldk + np.max([0, diverge.inflow]) - np.max([0, diverge.outflow])
            self.k = self.oldk + np.max([0, self.inflow]) - np.max([0, self.outflow])
            next_c.updated, self.updated, diverge.updated = True, True, True
                    
        else: # Normal cell
            if self.updated:
                return
            
            if len(self.cfrom) == 0:
                self.inflow = np.min([self.qmax, self.arr_rate, self.w * (self.kjam - self.oldk)]) * self.time_hour / self.length
                self.outflow = np.min([self.cto[0].qmax, self.oldk * self.vf, self.cto[0].w * (self.cto[0].kjam - self.cto[0].oldk)]) * self.time_hour / self.length
                
            elif len(self.cto) == 0:
                self.inflow = np.min([self.qmax, self.cfrom[0].oldk * self.vf, self.w * (self.kjam - self.oldk)]) * self.time_hour / self.length
                self.outflow = np.min([self.qmax, self.oldk * self.vf, self.dis_rate]) * self.time_hour / self.length
                
            else:
                self.inflow = np.min([self.qmax, self.cfrom[0].oldk * self.vf, self.w * (self.kjam - self.oldk)]) * self.time_hour / self.length
                self.outflow = np.min([self.cto[0].qmax, self.oldk * self.vf, self.cto[0].w * (self.cto[0].kjam - self.cto[0].oldk)]) * self.time_hour / self.length
            
            self.k = self.oldk + np.max([0, self.inflow]) - np.max([0, self.outflow])            
            self.updated = True

class interSection(object):
    idcase = {}
    count = 0
    def __init__(self, cell):
        if cell.__class__ != Cell:
            raise Exception("Input parameter must be a Cell class object that has more than one upstream or downstream cell.")
        self.id = interSection.count
        interSection.count += 1
        interSection.idcase[cell] = self
        
    def setProportion(self, p1, p2):
        self.p1 = p1 # probability from upstream normal cell or toward downstream normal cell 
        self.p2 = p2 # probability from upstream merge cell or toward downstream diverge cell 
    
    def getProportion(self):
        return (self.p1, self.p2)
    
class node(object):
    idcase = {}
    def __init__(self, nid, x, y, ntype):
        self.id = nid
        self.x = x
        self.y = y
        self.type = ntype
        node.idcase[nid] = self
        
    def getNodeFromID(nid):
        return node.idcase[nid]
    
    def getNodeByType(nodetype):
        need_return = []
        for key in node.idcase:
            if node.idcase[key].type == nodetype:
                need_return.append(node.idcase[key])
                
        return need_return
    
class link(object):
    idcase = {}
    def __init__(self, lid, fnode, tnode, speed, num_of_lanes, length):
        self.id = str(lid)
        self.source = str(fnode)
        self.sink = str(tnode)
        self.length = length
        self.speed = speed
        self.num_of_lanes = num_of_lanes
        link.idcase[str(lid)] = self
        
    def getLinkFromID(lid):
        return link.idcase[lid]
        
class trafficSignal(object):
    idcase = {}
    def __init__(self, tsid, corr_node):
        self.id = tsid
        self.node = corr_node
        self.phase = None
        self.turning = {}
        self.duration = 0
        trafficSignal.idcase[corr_node] = self
        
    def getPhaseByNode(corr_node):
        return trafficSignal.idcase[corr_node]

    def initPhase(self):
        corr_node = node.getNodeFromID(self.node)
        nodelist = []
        link_vector = {}
        if corr_node.type != 0:
            return 0

        for key in link.idcase:
            if link.idcase[key].sink != corr_node.id:
                continue

            nodelist.append(node.getNodeFromID(link.idcase[key].source))
            link_vector[key] = [node.getNodeFromID(link.idcase[key].source).x - corr_node.x, node.getNodeFromID(link.idcase[key].source).y - corr_node.y]

        for subkey in link_vector:
            if link_vector[subkey][1] >= link_vector[subkey][0]:
                if link_vector[subkey][1] >= -1 * link_vector[subkey][0]:
                    self.turning['southbound'] = link.getLinkFromID(subkey).source
                else:
                    self.turning['eastbound'] = link.getLinkFromID(subkey).source
            else:
                if link_vector[subkey][1] >= -1 * link_vector[subkey][0]:
                    self.turning['westbound'] = link.getLinkFromID(subkey).source
                else:
                    self.turning['northbound'] = link.getLinkFromID(subkey).source


    def setPhase(self, phase_num, max_phases):
        base_duration = 60
        left_turning_coef = 1/3
        sum_duration = 0

        self.phase = pd.DataFrame(columns=['id', 'phase', 'duration'])
        for i in range(max_phases):
            if (i + 1) % phase_num:
                duration = base_duration
            else:
                duration = base_duration * left_turning_coef

            if i < max_phases / phase_num and (i + 1) % phase_num:
                self.phase.loc[i] = [i, {self.turning.get("eastbound"):self.turning.get("westbound"), self.turning.get("westbound"):self.turning.get("eastbound")}, duration]
            elif i < max_phases / phase_num and not (i + 1) % phase_num:
                self.phase.loc[i] = [i, {self.turning.get("eastbound"):self.turning.get("southbound"), self.turning.get("westbound"):self.turning.get("northbound")}, duration]
            elif i >= max_phases / phase_num and (i + 1) % phase_num:
                self.phase.loc[i] = [i, {self.turning.get("southbound"):self.turning.get("northbound"), self.turning.get("northbound"):self.turning.get("southbound")}, duration]
            elif i >= max_phases / phase_num and not (i + 1) % phase_num:
                self.phase.loc[i] = [i, {self.turning.get("northbound"):self.turning.get("eastbound"), self.turning.get("southbound"):self.turning.get("westbound")}, duration]

            sum_duration += duration
            
        self.duration = sum_duration
        
    def getCurrentPhase(self, curr_time):
        curr_duration = curr_time % self.duration
        for i in range(len(self.phase)):
            curr_duration = curr_duration - self.phase.iloc[i]['duration']
            if curr_duration <= 0:
                return {'phase_timing':self.phase.iloc[i]['phase'], 'phase_number':self.phase.iloc[i]['id']}
    
def getCrossProduct(va, vb):
    return va[0]*vb[1] - va[1]*vb[0]

def getEuclideanDis(x1, x2, y1, y2):
    return np.sqrt(np.power(x2 - x1, 2) + np.power(y2 - y1, 2))

def readNetwork(filename):
    linkdf = pd.read_csv(filename, dtype={'id':np.str, 'from':np.str, 'to':np.str})
    nodedf = pd.read_csv("input_node.csv", dtype={'id':np.str})
    max_cell_length = 250 # meters
    junctions = []
    count = 0
    
    for ndidx in range(len(nodedf)):
        node(*list(nodedf.iloc[ndidx]))
    
    for index in range(len(linkdf)):
        link_length = getEuclideanDis(node.getNodeFromID(linkdf.iloc[index]['from']).x, \
                                      node.getNodeFromID(linkdf.iloc[index]['to']).x, \
                                      node.getNodeFromID(linkdf.iloc[index]['from']).y, \
                                      node.getNodeFromID(linkdf.iloc[index]['to']).y)
            
        link(*list(linkdf.iloc[index].drop("corresponding_link")), link_length)
            
        cell_length = link_length / (int(link_length / max_cell_length) + 1) / 1000
        cells = []
            
        for subidx in range(int(link_length / max_cell_length) + 1):
            link_level = 1
            cells.append(Cell('C'+str(subidx), str(int(linkdf.iloc[index]['id'])), 'A0', \
                 qmax=1800 / (int(int(link_level) / 3) + 1), kjam=120 / (int(int(link_level) / 3) + 1), \
                 vf=int(linkdf.iloc[index]['speed']), w=int(linkdf.iloc[index]['speed']) * 0.2, \
                 arr_rate=0, dis_rate=1800, length=cell_length))
                
        for subidx in range(len(cells)):
            if subidx < len(cells) - 1:
                cells[subidx].addConnection(cells[subidx + 1])
                
    for ndidx in range(len(nodedf)):
        try:
            if int(nodedf.iloc[ndidx]['id']) not in list(linkdf['from']):
                continue
            if nodedf.iloc[ndidx]['type'] == 0:
                junctions.append(trafficSignal(count, nodedf.iloc[ndidx]['id']))
                count += 1
        except:
            continue
            
    for elem in junctions:
        elem.initPhase()
        elem.setPhase(2, 4)
        
def readNetwork_NEXTA_format():
    newlinkdf = pd.read_csv("input_link_NEXTA.csv", dtype={'road_link':np.str, 'from_node_id':np.str, 'to_node_id':np.str})
    newnodedf = pd.read_csv("input_node_NEXTA.csv", dtype={'node_id':np.str})
    max_cell_length = 250 # meters
    
    linkdf = pd.DataFrame(columns=['id', 'from', 'to', 'speed', 'num_of_lanes', 'capacity'])
    nodedf = pd.DataFrame(columns=['id', 'x', 'y', 'type'])
    
    nodedf['id'] = newnodedf['node_id']
    nodedf['x'] = newnodedf['x_coord']
    nodedf['y'] = newnodedf['y_coord']
    nodedf['type'] = newnodedf['node_type']
    
    linkdf['id'] = newlinkdf['road_link_id']
    linkdf['from'] = newlinkdf['from_node_id']
    linkdf['to'] = newlinkdf['to_node_id']
    linkdf['speed'] = newlinkdf['free_speed']
    linkdf['num_of_lanes'] = newlinkdf['lanes']
    
    for ndidx in range(len(nodedf)):
        node(*list(nodedf.iloc[ndidx]))
    
    for index in range(len(linkdf)):
        link_length = newlinkdf.where(linkdf.iloc[index]['id'] == newlinkdf['road_link_id']).dropna(subset=['road_link_id']).iloc[0]['length']
        link_capacity = newlinkdf.where(linkdf.iloc[index]['id'] == newlinkdf['road_link_id']).dropna(subset=['road_link_id']).iloc[0]['capacity'] * linkdf.iloc[index]['num_of_lanes']
            
        link(*list(linkdf.iloc[index].drop("capacity")), link_length)
            
        cell_length = link_length / (int(link_length / max_cell_length) + 1) / 1000
        cells = []
            
        for subidx in range(int(link_length / max_cell_length) + 1):
            link_level = 1
            cells.append(Cell('C'+str(subidx), str(linkdf.iloc[index]['id']), 'A0', \
                 qmax=link_capacity, kjam=120 / (int(int(link_level) / 3) + 1), \
                 vf=int(linkdf.iloc[index]['speed']), w=int(linkdf.iloc[index]['speed']) * 0.2, \
                 arr_rate=0, dis_rate=link_capacity, length=cell_length))
                
        for subidx in range(len(cells)):
            if subidx < len(cells) - 1:
                cells[subidx].addConnection(cells[subidx + 1])
                
#    for key in link.idcase:
#        for subkey in link.idcase:
#            if link.idcase[key].source == link.idcase[subkey].sink:
#                Cell.getLastCell(subkey).addConnection(Cell.getFirstCell(key))
#            elif link.idcase[key].sink == link.idcase[subkey].source:
#                Cell.getLastCell(key).addConnection(Cell.getFirstCell(subkey))
                
    return 0
        
def quicklyCreateCells(number, linkid):
    cells = []
    for i in range(number):
        cells.append(Cell('C'+str(i), linkid, 'A0', qmax=5400, kjam=300, vf=80, w=16, arr_rate=1000, dis_rate=5400, length=0.2))
#        cells.append(Cell('C'+str(i), linkid, 'A0',arr_rate=1000, dis_rate=5400))
                
    for index in range(len(cells)):
        if index < len(cells) - 1:
            cells[index].addConnection(cells[index + 1])
                
    return cells

def timeDependentDemand(order=1, t=0, miu=1800, gamma=1, t0=0, t3=0, m=0.58):
    t = t / 600
    t0 = t0 / 600
    t3 = t3 / 600
    t2 = m * (t3 - t0) + t0
    if t < t0 or t > t3:
        return miu * 0.75
    if order == 1:
        return gamma * t + t0 + miu
    elif order == 2:
        return gamma * (t - t0)*(t2 - t) + miu
    elif order == 3:
#        tbar = t0 + (3*(t3 - t0)**2 - 4*(t2-t0)*(t3-t0)) / (4*(t3-t0) - 6*(t2-t0))
        tbar = t0 + (3 * (t3 - t0) - 4 * (t2 - t0)) / (4 - 6 * m)
        return miu + gamma * (t - t0)*(t - t2)*(t - tbar)
    else:
        raise Exception("Invaild input parameter! Order of time dependtent demand formula must be 1, 2 or 3")
        
        
def calibration():
    df = pd.read_csv("calibration.csv", dtype={'id':np.str})
    for key in Cell.idcase:
        cell = Cell.idcase[key]
        cell.qmax = df.where(df['id']==cell.linkid).dropna().iloc[0]['qmax']
        cell.kjam = df.where(df['id']==cell.linkid).dropna().iloc[0]['kjam']
        cell.vf = df.where(df['id']==cell.linkid).dropna().iloc[0]['vf']
        cell.w = df.where(df['id']==cell.linkid).dropna().iloc[0]['w']
        
def getLinkPerformance(outflow, traveltime):
    linklist = link.idcase.keys()
    df = pd.DataFrame(columns=['road_link_id', 'time_period', 'volume', 'travel_time'])
    newdf = df.copy()
    for linkid in linklist:
        cells = []
        for key in Cell.idcase:
            if Cell.idcase[key].linkid == linkid:
                cells.append(Cell.idcase[key].getCompleteAddress())
           
        df['volume'] = outflow.loc[cells].mean()
        df['travel_time'] = 1 / (traveltime.loc[cells].mean() / link.idcase[linkid].length / 1000 * 3600)
        df['time_period'] = outflow.loc[cells].mean().index
        temp_linkid_list = []
        for i in range(len(outflow.loc[cells].mean())):
            temp_linkid_list.append(linkid)
        df['road_link_id'] = temp_linkid_list
        newdf = newdf.append(df)
        
    return newdf
    
def simulation_Main(endtime):
    np.random.seed(123)
    newCellCase = {}
    dfindex = []
    
    Cell.getLastCell("0").addConnection(Cell.getFirstCell("1"))
    Cell.getLastCell("1").addConnection(Cell.getFirstCell("2"))
    miu_calibrated = {"0":1133.38, "1":1033.30, "2":849.02}
    
    Cell("C0", "M0", "A0", qmax=1200, kjam=120, vf=60, w=12, arr_rate=100)
    Cell("C0", "M1", "A0", qmax=1200, kjam=120, vf=60, w=12, arr_rate=300)
    Cell("C0", "M2", "A0", qmax=1200, kjam=120, vf=60, w=12, arr_rate=1200)
    Cell("C0", "D0", "A0", qmax=1200, kjam=120, vf=60, w=12, dis_rate=1200)
    Cell("C0", "D1", "A0", qmax=1200, kjam=120, vf=60, w=12, dis_rate=1200)
    Cell("C0", "D2", "A0", qmax=1200, kjam=120, vf=60, w=12, dis_rate=1200)
    Cell("C0", "D3", "A0", qmax=1200, kjam=120, vf=60, w=12, dis_rate=1200)
    
    Cell.getCell("A0.M0.C0").addConnection(Cell.getFirstCell("0").cto[0])
    Cell.getCell("A0.M1.C0").addConnection(Cell.getFirstCell("1").cto[0])
    Cell.getCell("A0.M2.C0").addConnection(Cell.getFirstCell("2").cto[0])
    
    Cell.getLastCell("0").cfrom[0].addConnection(Cell.getCell("A0.D0.C0"))
    Cell.getLastCell("1").cfrom[0].addConnection(Cell.getCell("A0.D1.C0"))
    Cell.getCell("A0.2.C4").addConnection(Cell.getCell("A0.D2.C0"))
    Cell.getLastCell("2").cfrom[0].addConnection(Cell.getCell("A0.D3.C0"))
    
    Cell.getCell("A0.2.C4").pk = 0.95
    Cell.getLastCell("2").cfrom[0].pk = 0.95
    
    for key in Cell.idcase:
        if Cell.idcase[key].cellid == 'C0' and 'M' not in Cell.idcase[key].linkid and 'D' not in Cell.idcase[key].linkid:
            if node.getNodeFromID(link.idcase[Cell.idcase[key].linkid].source).type == 6:
                Cell.idcase[key].arr_rate = Cell.idcase[key].qmax * 0.75
 
        newCellCase[Cell.idcase[key]] = Cell.idcase[key].cellid
    
#    newCellCase = dict(sorted(newCellCase.items(), key=lambda item:item[1]))
    for ncckey in newCellCase:
        dfindex.append(ncckey.getCompleteAddress())

    df = pd.DataFrame(index=dfindex)
    dfinflow = df.copy()
    dfoutflow = df.copy()
    
    t0_preset = 30
    t3_preset = 2580
                
    for t in range(endtime):
        density = []
        inflow = []
        outflow = []
        
        for cellelem in newCellCase:
            if cellelem.getCompleteAddress() == 'A0.0.C0':
                cellelem.arr_rate = timeDependentDemand(order=3, t = t, gamma=3.37, t0 = t0_preset, t3 = t3_preset, miu=993.79)
                
            if t >= t0_preset and t <= t3_preset and 'M' not in cellelem.linkid and 'D' not in cellelem.linkid and '2' in cellelem.linkid:
                cellelem.qmax = miu_calibrated[cellelem.linkid]
                
            cellelem.updateDensity()
            if cellelem.k <= 0:
                cellelem.k = 0
                
        for elem in newCellCase:
            density.append(elem.k)
            inflow.append(elem.inflow * elem.length / elem.time_hour)
            outflow.append(elem.outflow * elem.length / elem.time_hour)
            elem.updated = False
            
        df["t%i"%t] = density
        dfinflow["t%i"%t] = inflow
        dfoutflow["t%i"%t] = outflow
        
    return (df.sort_index(), dfinflow.sort_index(), dfoutflow.sort_index())
    
if __name__ == '__main__':
#    readNetwork_NEXTA_format()
    readNetwork("input_link.csv")
    calibration()
    output = simulation_Main(2600)
    output[0].to_csv("Density profile.csv")
    output[1].to_csv("Inflow profile.csv")
    output[2].to_csv("Outflow profile.csv")
    speed = (output[2] / output[0]).fillna(0)
    speed.to_csv("Speed profile.csv")
    lp = getLinkPerformance(output[2], speed)
    lp.to_csv("link_performance.csv")