#ifndef DIPOLECELL_MODELTISSUEPARAMETER_H
#define DIPOLECELL_MODELTISSUEPARAMETER_H

#include "CommonTypes.h"

namespace ModelParameter {
inline constexpr TissueParameter tissue_parameter = {
    .dt_full_step = 0.01,
    .dt_odeint = 0.001,
};
inline constexpr std::size_t num_pre_relaxation_steps = 500;
}
#endif // DIPOLECELL_MODELTISSUEPARAMETER_H
