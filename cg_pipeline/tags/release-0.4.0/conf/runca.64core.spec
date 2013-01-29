cnsErrorRate = 0.10
ovlErrorRate = 0.10
 
overlapper = ovl
unitigger = bogart
utgBubblePopping = 1
 
merSize = 10
 
merylMemory = 12800
merylThreads = 16
 
ovlStoreMemory = 16384
 
# grid info
useGrid = 0
scriptOnGrid = 0
frgCorrOnGrid = 1
ovlCorrOnGrid = 1
 
sge = -A assembly
sgeScript = -pe threads 16
sgeConsensus = -pe threads 1
sgeOverlap = -pe threads 2
sgeFragmentCorrection = -pe threads 2
sgeOverlapCorrection = -pe threads 1
 
#ovlMemory=8GB --hashload 0.7
ovlHashBits=25
ovlThreads = 20
ovlHashBlockLength = 30000000
ovlRefBlockSize =  5000000
 
# for mer overlapper
merCompression = 1
merOverlapperSeedBatchSize = 500000
merOverlapperExtendBatchSize = 250000
 
frgCorrThreads = 20
frgCorrBatchSize = 100000
 
ovlCorrBatchSize = 100000
 
# non-Grid settings, if you set useGrid to 0 above these will be used
merylMemory = 12800
merylThreads = 20
 
ovlStoreMemory = 16384
 
ovlConcurrency = 40
 
merOverlapperThreads = 20
merOverlapperSeedConcurrency = 40
merOverlapperExtendConcurrency = 40
 
frgCorrConcurrency = 40
 
ovlCorrConcurrency = 40
cnsConcurrency = 40
 
doToggle=0
toggleNumInstances = 0
toggleUnitigLength = 2000
 
 
 
doOverlapBasedTrimming = 1
doExtendClearRanges = 2


