#pragma once
// #if defined(_WIN32)
// #include <windows.h>

#include <vector>
#include <math.h>
#include <cstring>
#include "osutils.h"

namespace aadc {

	class JitRuntime {
		// pointer and number of bytes allocated
		std::vector<std::pair<void*, size_t>> _data;
	public:
		JitRuntime() {
		}
		~JitRuntime() {
			// release all allocated blocks
			for (auto pair : _data) {
				osutils::releaseMem(pair.first, pair.second);
			}
		}

		template<typename Func>
		void Add(Func* dst, const std::vector<unsigned char>& code) {
			auto meminfo   = osutils::getMemInfo();
			auto page_size = meminfo.pageSize;

			if (page_size < code.size()) {
				auto pages = (uint32_t)::floor(code.size() / page_size)+1;
				page_size  = pages * page_size;
			}

			auto const buffer = osutils::allocMem(page_size);
			if (buffer == nullptr)
				throw std::runtime_error("Not able to allocate memory");

			_data.push_back(std::make_pair(buffer, page_size));

			// copy the machine code into that memory:
			std::memcpy(buffer, code.data(), code.size());

			// mark memory as safe to execute
			osutils::protectMem(buffer, code.size());

			*dst = reinterpret_cast<Func>(buffer);
		}
	};
}

// #endif
