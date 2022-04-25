#pragma once
#include "Bone.h"
#include <vector>
//a skeleton is just a list of Bones. Each bone knows its parent as an index into this list.
using Skeleton = std::vector<Bone>;

