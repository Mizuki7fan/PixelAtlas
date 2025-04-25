#include "process.h"
#include "naive_pixelator.h"
#include <frame/global_defs.h>
#include <frame/io.h>
#include <frame/metric.h>

namespace fs = std::filesystem;
void MainProcess() { fs::path instance_path = frm::global::InstancePath(); }