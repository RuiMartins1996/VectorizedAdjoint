#pragma once

#include <aadc/api.h>
#include <memory>

namespace aadc {

	struct MemInfo {
		size_t pageSize;		// Virtual memory page size.
		bool   cpu_AVX512;		// cpu supports AVX512
	};
};

extern "C" AADC_API const aadc::MemInfo* CAAD_OSUtilsGetMemInfo();
extern "C" AADC_API void* CAAD_AllocMem(size_t size);
extern "C" AADC_API void CAAD_ProtectMem(void* p, size_t size);
extern "C" AADC_API void CAAD_ReleaseMem(void* p, size_t size);

namespace aadc {

	struct osutils {

		static const MemInfo& getMemInfo() {
			return *(CAAD_OSUtilsGetMemInfo());
		}

		static void* allocMem(size_t size) {
			return CAAD_AllocMem(size);
		}
		static void  protectMem(void* p, size_t size) {
			CAAD_ProtectMem(p,size);
		}
		static void  releaseMem(void* p, size_t size) {
			CAAD_ReleaseMem(p, size);
		}
	};

}