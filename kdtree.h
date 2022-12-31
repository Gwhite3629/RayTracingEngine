#ifndef _KDTREE_H_
#define _KDTREE_H_

#include <algorithm>
#include <cmath>
#include <vector>

template <typename T>
struct data {
    T containter;
    std::vector<int> p;
};

template <typename T>
class KDTree
{
public:
    explicit KDTree(const std::vector<data<T>> &points)
        : num_cells_(std::cbrt(points.size())),
          cells_(num_cells_,
                 std::vector<std::vector<std::vector<data<T>>>>(
                     num_cells_, std::vector<std::vector<data<T>>>(num_cells_)))
    {
        // Determine the minimum and maximum values along each dimension.
        std::vector<int> min(3, std::numeric_limits<int>::max());
        std::vector<int> max(3, std::numeric_limits<int>::min());
        for (const auto &point : points)
        {
            for (int i = 0; i < 3; i++)
            {
                min[i] = std::min(min[i], point.p[i]);
                max[i] = std::max(max[i], point.p[i]);
            }
        }
        min_ = min;
        max_ = max;

        // Calculate the cell size along each dimension.
        std::vector<int> cell_size(3);
        for (int i = 0; i < 3; i++)
        {
            cell_size[i] = (max[i] - min[i]) / num_cells_;
        }

        cell_size_ = cell_size;

        // Add the points to the appropriate cells.
        for (const auto &point : points)
        {
            // Calculate the cell indices for the point.
            std::vector<int> cell_indices(3);
            for (int i = 0; i < 3; i++)
            {
                cell_indices[i] = std::min((point.p[i] - min[i]) / cell_size[i], num_cells_ - 1);
            }

            // Add the point to the appropriate cell.
            cells_[cell_indices[0]][cell_indices[1]][cell_indices[2]].push_back(point);
        }
    }

    void insert(const data<T> &point)
    {
        // Calculate the cell indices for the point.
        std::vector<int> cell_indices(3);
        for (int i = 0; i < 3; i++)
        {
            cell_indices[i] = std::min((point.p[i] - min_[i]) / cell_size_[i], num_cells_ - 1);
        }

        // Add the point to the appropriate cell.
        cells_[cell_indices[0]][cell_indices[1]][cell_indices[2]].push_back(point);
    }

    data<T> nearestNeighbor(const std::vector<int> &point) const
    {
        // Calculate the cell indices for the point.
        std::vector<int> cell_indices(3);
        for (int i = 0; i < 3; i++)
        {
            cell_indices[i] = std::min((point[i] - min_[i]) / cell_size_[i], num_cells_ - 1);
        }

        // Start with the cell containing the point as the best candidate.
        data<T> nearest;
        nearest.p = point;
        int best_distance = std::numeric_limits<int>::max();

        // Search the cells in a spiral pattern starting from the cell containing the point.
        int dx[] = {0, 1, 0, -1};
        int dy[] = {1, 0, -1, 0};
        int dz[] = {0, 0, 0, 0};
        int x = cell_indices[0];
        int y = cell_indices[1];
        int z = cell_indices[2];
        int dir = 0;
        for (int i = 0; i < num_cells_; i++)
        {
            for (int j = 0; j < i + 1; j++)
            {
                // Check the points in the current cell.
                for (const auto &candidate : cells_[x][y][z])
                {
                    int distance = squaredDistance(point, candidate.p);
                    if (distance < best_distance)
                    {
                        nearest = candidate;
                        best_distance = distance;
                    }
                }

                // Move to the next cell in the spiral pattern.
                x += dx[dir];
                y += dy[dir];
                z += dz[dir];
                if (dir == 0 && y == cell_indices[1] + i + 1)
                {
                    dir = 1;
                }
                else if (dir == 1 && x == cell_indices[0] - i - 1)
                {
                    dir = 2;
                }
                else if (dir == 2 && y == cell_indices[1] - i - 1)
                {
                    dir = 3;
                }
                else if (dir == 3 && x == cell_indices[0] + i + 1)
                {
                    dir = 0;
                }
            }
        }

        return nearest;
    }

    std::vector<data<T>> rangeSearch(const std::vector<int> &lower,
                                            const std::vector<int> &upper) const
    {
        std::vector<data<T>> points_in_range;

        // Calculate the cell indices for the lower and upper bounds of the range.
        std::vector<int> lower_indices(3);
        std::vector<int> upper_indices(3);
        for (int i = 0; i < 3; i++)
        {
            lower_indices[i] = std::max(0, (lower[i] - min_[i]) / cell_size_[i]);
            upper_indices[i] = std::min(num_cells_ - 1, (upper[i] - min_[i]) / cell_size_[i]);
        }

        // Search the cells within the range.
        for (int x = lower_indices[0]; x <= upper_indices[0]; x++)
        {
            for (int y = lower_indices[1]; y <= upper_indices[1]; y++)
            {
                for (int z = lower_indices[2]; z <= upper_indices[2]; z++)
                {
                    // Check the points in the current cell.
                    for (const auto &point : cells_[x][y][z])
                    {
                        if (point.p[0] >= lower[0] && point.p[0] <= upper[0] &&
                            point.p[1] >= lower[1] && point.p[1] <= upper[1] &&
                            point.p[2] >= lower[2] && point.p[2] <= upper[2])
                        {
                            points_in_range.push_back(point);
                        }
                    }
                }
            }
        }

        return points_in_range;
    }

    std::vector<data<T>> kNearestNeighbors(const std::vector<int> &point,
                                                  int k) const
    {
        std::vector<std::pair<int, data<T>>> distances;

        // Calculate the cell indices for the point.
        std::vector<int> cell_indices(3);
        for (int i = 0; i < 3; i++)
        {
            cell_indices[i] = std::min((point[i] - min_[i]) / cell_size_[i], num_cells_ - 1);
        }

        // Search the cells in a spiral pattern starting from the cell containing the point.
        int dx[] = {0, 1, 0, -1};
        int dy[] = {1, 0, -1, 0};
        int dz[] = {0, 0, 0, 0};
        int x = cell_indices[0];
        int y = cell_indices[1];
        int z = cell_indices[2];
        int dir = 0;
        for (int i = 0; i < num_cells_; i++)
        {
            for (int j = 0; j < i + 1; j++)
            {
                // Check the points in the current cell.
                for (const auto &candidate : cells_[x][y][z])
                {
                    int distance = squaredDistance(point, candidate.p);
                    distances.emplace_back(distance, candidate.p);
                }

                // Move to the next cell in the spiral pattern.
                x += dx[dir];
                y += dy[dir];
                z += dz[dir];
                if (dir == 0 && y == cell_indices[1] + i + 1)
                {
                    dir = 1;
                }
                else if (dir == 1 && x == cell_indices[0] - i - 1)
                {
                    dir = 2;
                }
                else if (dir == 2 && y == cell_indices[1] - i - 1)
                {
                    dir = 3;
                }
                else if (dir == 3 && x == cell_indices[0] + i + 1)
                {
                    dir = 0;
                }
            }
        }
        // Sort the distances and return the k nearest neighbors.
        std::sort(distances.begin(), distances.end());
        std::vector<data<T>> nearest_neighbors;
        for (int i = 0; i < k && i < distances.size(); i++)
        {
            nearest_neighbors.push_back(distances[i].p.second);
        }

        return nearest_neighbors;
    }

    void clear()
    {
        for (int x = 0; x < num_cells_; x++)
        {
            for (int y = 0; y < num_cells_; y++)
            {
                for (int z = 0; z < num_cells_; z++)
                {
                    cells_[x][y][z].p.clear();
                }
            }
        }
    }

    int size() const
    {
        int size = 0;
        for (int x = 0; x < num_cells_; x++)
        {
            for (int y = 0; y < num_cells_; y++)
            {
                for (int z = 0; z < num_cells_; z++)
                {
                    size += cells_[x][y][z].size();
                }
            }
        }
        return size;
    }

    bool isEmpty() const
    {
        for (int x = 0; x < num_cells_; x++)
        {
            for (int y = 0; y < num_cells_; y++)
            {
                for (int z = 0; z < num_cells_; z++)
                {
                    if (!cells_[x][y][z].p.empty())
                    {
                        return false;
                    }
                }
            }
        }
        return true;
    }

    

private:
    // The number of cells in the divided space.
    int num_cells_;

    // The cells containing the points.
    std::vector<std::vector<std::vector<std::vector<data<T>>>>> cells_;

    std::vector<int> min_;
    std::vector<int> max_;
    std::vector<int> cell_size_;

    static int squaredDistance(const std::vector<int> &p1,
                             const std::vector<int> &p2)
    {
        int squared_distance = 0;
        for (int i = 0; i < 3; ++i)
        {
            squared_distance += (p1[i] - p2[i]) * (p1[i] - p2[i]);
        }
        return squared_distance;
    }
};

#endif // _KDTREE_H_