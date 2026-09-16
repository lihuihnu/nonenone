/**
 * @file geometry.hpp
 * @brief CpGrid 几何量计算与向量辅助。
 */
#pragma once

#include <cpgrid/point.hpp>

#include <vector>

namespace MPMC
{

class Node;

/**
 * @brief 计算三维叉积。
 */
[[nodiscard]] Point cross(
    const Point &a,
    const Point &b) noexcept;

/**
 * @brief 计算三角形面积。
 */
[[nodiscard]] double triangleArea(
    const Point &a,
    const Point &b,
    const Point &c) noexcept;

/**
 * @brief 计算四面体体积绝对值。
 */
[[nodiscard]] double tetrahedronVolume(
    const Point &a,
    const Point &b,
    const Point &c,
    const Point &d) noexcept;

/**
 * @brief 计算节点集合的算术平均点。
 */
[[nodiscard]] Point averagePoint(
    const std::vector<Node *> &nodes);

/**
 * @brief 一次多边形遍历得到面积、质心和单位法向量。
 *
 * 该结果用于 CpGrid setup 的融合几何内核，避免对同一 face 为 area /
 * centroid / normal 重复三次参考点和边循环。
 */
struct PolygonGeometry final
{
    double area{0.0};
    Point centroid{};
    Point unitNormal{};
};

/**
 * @brief 一次遍历计算平面多边形面积、面积加权中心和单位法向量。
 *
 * 数值公式与 polygonArea()/polygonCentroid()/polygonUnitNormal() 保持一致，
 * 只是共享 reference point 和 edge traversal。
 */
[[nodiscard]] PolygonGeometry polygonGeometry(
    const std::vector<Node *> &nodes,
    const Point &referencePoint);

/**
 * @brief 计算平面多边形面积。
 *
 * 使用由 referencePoint 到相邻顶点组成的三角形扇区求和。
 */
[[nodiscard]] double polygonArea(
    const std::vector<Node *> &nodes,
    const Point &referencePoint);

/**
 * @brief 计算平面多边形的面积加权中心。
 */
[[nodiscard]] Point polygonCentroid(
    const std::vector<Node *> &nodes,
    const Point &referencePoint);

/**
 * @brief 计算多边形单位法向量。
 *
 * 使用三角扇区面积向量累加，方向由节点顺序决定。
 */
[[nodiscard]] Point polygonUnitNormal(
    const std::vector<Node *> &nodes,
    const Point &referencePoint);

} // namespace MPMC
