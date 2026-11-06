// Export a serial Cartesian AMReX plotfile without unit conversion.
#include <AMReX.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_MultiFab.H>
#include <AMReX_BoxIterator.H>
#include <AMReX_ParallelDescriptor.H>
#include <fstream>
#include <iomanip>
#include <map>
#include <tuple>
#include <vector>

int main(int argc, char* argv[])
{
    if (argc != 3 && argc != 4) { return 2; }
    const bool composite = argc == 4 && std::string(argv[3]) == "--composite";
    const bool all_levels = argc == 4 && std::string(argv[3]) == "--levels";
    const bool hierarchy = composite || all_levels;
    if (argc == 4 && !hierarchy) { return 2; }
    const std::string input = argv[1], output = argv[2];
    amrex::Initialize(argc, argv, false);
    {
        AMREX_ALWAYS_ASSERT(amrex::ParallelDescriptor::NProcs() == 1);
        amrex::PlotFileData plot(input);
        AMREX_ALWAYS_ASSERT(plot.spaceDim() == 3 && plot.coordSys() == 0);
        AMREX_ALWAYS_ASSERT(hierarchy || plot.finestLevel() == 0);
        const auto& names = plot.varNames();
        const auto lo = plot.probLo();
        std::ofstream stream(output);
        AMREX_ALWAYS_ASSERT(stream.good());
        if (hierarchy) stream << "level,dx,dy,dz,";
        stream << "time,x,y,z";
        for (const auto& name : names) stream << ',' << name;
        stream << '\n' << std::setprecision(17);
        for (int level = 0; level <= plot.finestLevel(); ++level) {
        auto field = plot.get(level);
        const auto dx = plot.cellSize(level);
        const auto domain = plot.probDomain(level);
        amrex::BoxArray covered;
        if (level < plot.finestLevel()) {
            covered = plot.boxArray(level + 1);
            covered.coarsen(plot.refRatioVect(level));
        }
        std::map<std::tuple<int,int,int>, std::vector<double>> rows;
        for (amrex::MFIter mfi(field); mfi.isValid(); ++mfi) {
            const auto values = field.const_array(mfi);
            for (amrex::BoxIterator cell(mfi.validbox()); cell.ok(); ++cell) {
                const auto index = cell();
                if (composite && covered.contains(index)) continue;
                std::vector<double> row;
                for (int n=0; n<field.nComp(); ++n) row.push_back(values(index[0],index[1],index[2],n));
                auto inserted = rows.emplace(std::make_tuple(index[2],index[1],index[0]), std::move(row));
                AMREX_ALWAYS_ASSERT(inserted.second);
            }
        }
        if (!hierarchy) AMREX_ALWAYS_ASSERT(rows.size() == static_cast<std::size_t>(domain.numPts()));
        if (all_levels) AMREX_ALWAYS_ASSERT(rows.size() == static_cast<std::size_t>(plot.boxArray(level).numPts()));
        for (const auto& [index,row] : rows) {
            const auto [k,j,i] = index;
            if (hierarchy) stream << level << ',' << dx[0] << ',' << dx[1] << ',' << dx[2] << ',';
            stream << plot.time() << ',' << lo[0]+(i-domain.smallEnd(0)+.5)*dx[0]
                   << ',' << lo[1]+(j-domain.smallEnd(1)+.5)*dx[1]
                   << ',' << lo[2]+(k-domain.smallEnd(2)+.5)*dx[2];
            for (double value : row) stream << ',' << value;
            stream << '\n';
        }
        }
        stream.close();
        AMREX_ALWAYS_ASSERT(stream.good());
    }
    amrex::Finalize();
}
