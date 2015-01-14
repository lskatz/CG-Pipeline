cnsErrorRate = 0.10
ovlErrorRate = 0.10
 
overlapper = ovl
unitigger = bogart
utgBubblePopping = 1
 
merSize = 14
 
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
ovlHashBits=23
ovlHashBlockLength= 30000000
ovlThreads = 2
ovlHashBlockLength = 20000000
ovlRefBlockSize =  5000000
 
# for mer overlapper
merCompression = 1
merOverlapperSeedBatchSize = 500000
merOverlapperExtendBatchSize = 250000
 
frgCorrThreads = 2
frgCorrBatchSize = 100000
 
ovlCorrBatchSize = 100000
 
# non-Grid settings, if you set useGrid to 0 above these will be used
merylMemory = 128000
merylThreads = 4
 
ovlStoreMemory = 8192
 
ovlConcurrency = 2
 
merOverlapperThreads = 3
merOverlapperSeedConcurrency = 2
merOverlapperExtendConcurrency = 2
 
frgCorrConcurrency = 2
 
ovlCorrConcurrency = 2
cnsConcurrency = 2
 
doToggle=0
toggleNumInstances = 0
toggleUnitigLength = 2000
 
doOverlapBasedTrimming = 1
doExtendClearRanges = 2
