import pandas as pd
import numpy as np

import copy
import math


class NeighborSearcher:
    def __init__(self, points, bounding, spacing=0.25):
        self.PointsNum = len(points)
        self.Points = points

        self.BoundingBox = bounding
        self.Bound = [(self.BoundingBox[1] - self.BoundingBox[0]),
                      (self.BoundingBox[3] - self.BoundingBox[2]),
                      (self.BoundingBox[5] - self.BoundingBox[4])]
        self.CellSize = spacing
        self.XYZCellNum = [math.ceil(self.Bound[0] / self.CellSize),
                           math.ceil(self.Bound[1] / self.CellSize),
                           math.ceil(self.Bound[2] / self.CellSize)]

        self.CellNum = self.XYZCellNum[0] * self.XYZCellNum[1] * self.XYZCellNum[2]

        self.HashList = [0] * self.PointsNum
        self.IndexList = [0] * self.PointsNum
        self.StartList = {}
        self.EndList = {}
        # self.numParInCell = [0] * self.CellNum
        # self.neighbor = [[] for i in range(self.numPar)]
        # self.numNeighbor = [0] * self.numPar
        # self.cell2Par = [[] for i in range(self.CellNum)]
        # self.numCell2Par = [0] * self.numPar
        self.init_hash_grid()

    def get_xyz_id(self, pos):
        return ((pos[0] - self.BoundingBox[0]) // self.CellSize,
                (pos[1] - self.BoundingBox[2]) // self.CellSize,
                (pos[2] - self.BoundingBox[4]) // self.CellSize)

    def get_cell_id(self, xyz_id):
        if (xyz_id[0] < 0 or xyz_id[0] >= self.XYZCellNum[0] or
                xyz_id[1] < 0 or xyz_id[1] >= self.XYZCellNum[1] or
                xyz_id[2] < 0 or xyz_id[2] >= self.XYZCellNum[2]):
            return -1
        return int(xyz_id[2] * self.XYZCellNum[0] * self.XYZCellNum[1] +
                   xyz_id[1] * self.XYZCellNum[0] + xyz_id[0])

    def init_hash_grid(self):
        # update hash table: cell2Par
        for i in range(self.PointsNum):
            cell_id = self.get_cell_id(self.get_xyz_id(self.Points[i]))  # get cell ID from position
            self.HashList[i] = cell_id
            self.IndexList[i] = i
        self.IndexList.sort(key=lambda x: self.HashList[x])
        temp = copy.deepcopy(self.HashList)
        for i in range(self.PointsNum):
            self.HashList[i] = temp[self.IndexList[i]]
        count = 0
        previous = -1
        for index in range(self.PointsNum):
            hash_id = self.HashList[index]
            if hash_id < 0:
                continue
            if hash_id != previous:
                self.StartList[hash_id] = count
                previous = hash_id
            count += 1
            self.EndList[hash_id] = count

    def get_in_cell_list(self, hash_id):
        pList = []
        if self.StartList.get(hash_id) is not None and self.EndList.get(hash_id) is not None:
            start_index = self.StartList[hash_id]
            end_index = self.EndList[hash_id]
        else:
            return pList
        for index in range(start_index, end_index):
            pList.append(self.IndexList[index])
        return pList
