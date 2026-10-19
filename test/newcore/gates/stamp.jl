# The commit a gate table is measured on, read when the script starts: the
# short hash of HEAD, and whether the working tree differed from it.
const COMMIT_STAMP = let repo = joinpath(@__DIR__, "..", "..", "..")
    head = readchomp(`git -C $repo rev-parse --short HEAD`)
    clean = isempty(readchomp(`git -C $repo status --porcelain`))
    clean ? "commit $(head)" : "commit $(head) with uncommitted changes"
end
