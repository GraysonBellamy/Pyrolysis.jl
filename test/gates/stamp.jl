# The commit a measurement table is measured on, read when the script starts:
# the short hash of HEAD and of its tree (the tree survives re-signing and
# rebasing, the commit hash does not), and whether the working tree differed
# from it outside the tables the measurement scripts write.
const COMMIT_STAMP = let repo = joinpath(@__DIR__, "..", "..")
    head = readchomp(`git -C $repo rev-parse --short HEAD`)
    tree = readchomp(`git -C $repo rev-parse --short HEAD^\{tree\}`)
    tables = ":!test/gates/measurements"
    clean = isempty(readchomp(`git -C $repo status --porcelain -- . $tables`))
    stamp = "commit $(head) (tree $(tree))"
    clean ? stamp : "$(stamp) with uncommitted changes"
end
