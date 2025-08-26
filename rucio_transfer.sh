rucio add-rule \
    cms:/SingleMuon/Run2016B-ZMu-21Feb2020_ver2_UL2016_HIPM-v1/RAW-RECO \
    1 \
    'rse_type=DISK&cms_type=real\tier=3\tier=0' \
    --lifetime 2592000 \
    --grouping 'ALL' \
    --activity "User AutoApprove" \
    --ask-approval \
    --comment "Details for use, ticket reference if any"
