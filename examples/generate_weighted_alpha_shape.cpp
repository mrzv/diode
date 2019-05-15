#include <vector>
#include <array>
#include <string>
#include <fstream>
#include <iostream>

#include <diode/diode.h>

struct Points
{
    float   operator()(unsigned i, unsigned j) const    { return points[i][j]; }
    size_t  size() const                                { return points.size(); }

    std::vector<std::array<float, 4>>   points;
};

struct OutputSimplex
{
    template<unsigned long D>
    void        operator()(const std::array<unsigned, D>& vertices, float alpha) const
    {
        for (unsigned i = 0; i < D; ++i)
            std::cout << vertices[i] << ' ';
        std::cout << alpha << std::endl;
    }
};

Points read_points(std::string fn)
{
    Points result;

    std::ifstream in(fn);
    while (in)
    {
        std::string line;
        std::getline(in, line);
        if (line[0] == '#') continue;

        std::istringstream iss(line);
        float x,y,z,w;
        iss >> x >> y >> z >> w;
        std::array<float, 4> p {x,y,z,w};
        result.points.emplace_back(p);
    }

    return result;
}

int main(int argc, char** argv)
{
    if (argc < 2)
    {
        std::cout << "Usage: " << argv[0] << " points" << std::endl;
        return 1;
    }

    // read points
    Points points = read_points(argv[1]);

    // generate the alpha shape
    diode::AlphaShapes<>::fill_weighted_alpha_shapes(points, OutputSimplex());
}

