# Release runbook

Step-by-step checklist for cutting a MirDeep-P3 release.  Written after the
3.1.6c release, so it records the order that actually works plus the traps that
cost time.  Everything here was executed and verified for 3.1.6c.

Throughout, `<VER>` is the bare version (e.g. `3.1.6c`) and `<REL>` is the
source tag (e.g. `v3.1.6c`).

---

## 0. What a release consists of

| Artifact | Produced by | Where it lands |
|---|---|---|
| Source tag `<REL>` | you | GitHub |
| Source tarball | GitHub, on demand | `github.com/.../archive/refs/tags/<REL>.tar.gz` |
| Docker image | `docker-release.yml`, on `<REL>-full` | Docker Hub `merc3dez/mirdeep-p3:<VER>-full`, GHCR `ghcr.io/yangxz-lab/mirdeep-p3:<VER>-full` |
| Offline image tarball | same workflow (`docker save \| gzip`) | GitHub Release `<REL>-full`, asset `mirdeep-p3-<VER>-full.tar.gz` |
| Aliyun mirror image | **manual** | `crpi-rs803yb7s70369gn.cn-beijing.personal.cr.aliyuncs.com/merc3dez/mirdeep-p3:<VER>-full` |
| conda package, personal channel | **manual** (`build_full_<VER>.sh`) | `anaconda.org/jaguares/mirdeep-p3` |
| conda package, bioconda | **manual** (PR to `bioconda-recipes`) | `anaconda.org/bioconda/mirdeep-p3` |

---

## 1. Pre-flight: make `main` green

```bash
R=/nfs/users/zhaojw/software/github/git/mirdeep-p3
cd $R && git status --short          # expect clean
git fetch origin && git status -sb   # expect no divergence
```

- CI runs on **Python 3.9**.  There is a ready local env (`test_mirdp3`) with all
  the external tools, so the CI path can be replayed locally:
  ```bash
  PY39=/nfs/users/zhaojw/miniconda3/envs/test_mirdp3/bin/python
  $PY39 mirdeep-p3 --version && $PY39 mirdeep-p3 -h
  $PY39 -c "import sys; sys.path.insert(0,'src'); import commands.identification, commands.annotation, commands.analysis"
  ```
- Push, then **wait for CI before tagging** (§9 shows how to query Actions
  without `gh`, which is not installed).

## 2. Bump the version — everywhere

`src/version.py` is the single source of truth, but the string is quoted in
several places.  Miss one and the release is inconsistent.

| File | What to change |
|---|---|
| `src/version.py` | `__version__` |
| `conda.recipe/meta.yaml` | `package.version` + `source.url` (leave `sha256` as a placeholder, §3 fills it) |
| `README.md` | Docker image names / download URLs / pin advice |
| `Dockerfile` | header comment example |
| `CHANGELOG.md` | add the new section |
| `CONTRIBUTING.md`, `test.sh`, `.github/workflows/*` | only if they quote a version.  The data-index tag comes from `data-index.env` and must **not** be touched |

Note there is **no version constant in `docker-release.yml`** — it derives the
image version from the tag (`VER="${GITHUB_REF_NAME#v}"`).

Verify nothing stale is left:
```bash
git grep -nE "3\.1\.[0-9][a-z]?" -- . | grep -v '^data/' | grep -v ko00001 | grep -v CHANGELOG.md
# only hits should be the data-index release URL, which is version-independent by design
```
Render the recipe to catch YAML/Jinja breakage:
```bash
export PATH=/nfs/users/zhaojw/miniconda3/bin:$PATH
conda render conda.recipe >/dev/null && echo "recipe renders"
```

## 3. Source tag, then the conda sha256

```bash
cd $R
git tag -a <REL> -m "MirDeep-P3 <VER>"        # annotated, matches history
git push origin <REL>

# the repo is public -> no token needed
curl -sL -o /tmp/<REL>.tar.gz \
  https://github.com/YangXZ-lab/mirdeep-p3/archive/refs/tags/<REL>.tar.gz
sha256sum /tmp/<REL>.tar.gz                    # -> paste into conda.recipe/meta.yaml
```
The digest is reproducible: `github.com/archive/...` and
`codeload.github.com/tar.gz/...` return the same file.  The tarball's top-level
directory is `mirdeep-p3-<VER>/`.

Commit the filled-in `sha256` and push.  (The recipe *inside* the tarball still
holds the placeholder — harmless: conda-build reads the working copy, and
bioconda carries the hash in its own recipe file.)

## 4. `-full` tag → image + GitHub Release

```bash
git tag -a <REL>-full -m "MirDeep-P3 <VER> full image"
git push origin <REL>-full     # triggers docker-release.yml
```

- Trigger is `tags: ["v*-full", "mirdeep-p3-v*"]`.  The second pattern is legacy
  (§8) — do **not** invent tags starting with `mirdeep-p3-v`.
- Convention: `<REL>-full` points at a **later** commit than `<REL>`
  (e.g. `v3.1.5a` = `c5025b3` while `v3.1.5a-full` = `58b8cb7`).  Tag
  `<REL>-full` at the tip of `main` after the sha256 commit.  Do not re-tag.
- The workflow downloads `data-index.tar.gz` from the Release named in
  `data-index.env`, builds, pushes to both registries, `docker save | gzip`s the
  image, and creates the GitHub Release with that tarball attached.

Confirm what was actually pushed.  The CI log is the only reliable source, since
`hub.docker.com` is unreachable from the cluster:
```bash
curl -sL -H "Authorization: token $CRED" \
  ".../actions/runs/$RID/logs" -o /tmp/l.zip
unzip -p /tmp/l.zip '*Build and push.txt' | grep -E "pushing manifest for|exporting config"
# expect two lines with the same digest: docker.io/... and ghcr.io/...
```

## 5. Offline image tarball (for Aliyun)

`compute1` cannot reach Docker Hub directly (`registry-1.docker.io` and
`auth.docker.io` time out).  `daemon.json` has domestic mirrors, but a 5.4 GB
pull takes >25 min and will not survive a foreground SSH timeout.

**Use the GitHub Release asset instead** — it *is* the `docker save | gzip`
output from CI, is only ~1.9 GB, and GitHub is directly reachable:
```bash
OUT=/nfs/users/zhaojw/Genburten/release
mkdir -p $OUT && cd $OUT
curl -fL --retry 3 -o mirdeep-p3-<VER>-full.tar.gz \
  https://github.com/YangXZ-lab/mirdeep-p3/releases/download/<REL>-full/mirdeep-p3-<VER>-full.tar.gz
sha256sum mirdeep-p3-<VER>-full.tar.gz
```

Verify it really is the published image before handing it over:
```bash
tar tzf mirdeep-p3-<VER>-full.tar.gz | grep -E "manifest.json|index.json"
tar xzf mirdeep-p3-<VER>-full.tar.gz manifest.json -O   # RepoTags must match
# and match the config blob name against the "exporting config sha256:..." line in the CI log
docker load -i mirdeep-p3-<VER>-full.tar.gz
docker run --rm merc3dez/mirdeep-p3:<VER>-full --version
```
> The image entrypoint is `/opt/mirdeep-p3/mirdeep-p3`, so it is
> `docker run <img> --version`, **not** `docker run <img> mirdeep-p3 --version`.

Upload (credentials are already stored on compute1 — no `docker login` needed):
```bash
docker tag merc3dez/mirdeep-p3:<VER>-full \
  crpi-rs803yb7s70369gn.cn-beijing.personal.cr.aliyuncs.com/merc3dez/mirdeep-p3:<VER>-full
docker push crpi-rs803yb7s70369gn.cn-beijing.personal.cr.aliyuncs.com/merc3dez/mirdeep-p3:<VER>-full
```
> The login username is the **Alibaba Cloud account name**, not the namespace
> `merc3dez`.  The ACR password is a self-set "fixed password" that the console
> never displays — if lost, it can only be reset (ACR console -> personal
> instance -> *Access Credentials*).

## 6. conda package for the personal channel (`jaguares`)

`conda.recipe/build.sh` copies `$SRC_DIR/data`, but `data/index` is `.gitignore`d
and therefore absent from the tag tarball.  A plain `conda build conda.recipe`
produces a package whose `identification` step cannot run.  Build it from a
**local path source** with the index placed next to it:

```bash
# Genburten/temp/build_full_<VER>.sh does exactly this, with a trap that always
# restores the URL recipe and removes the copied index:
#   1. back up the URL-source meta.yaml
#   2. move previously built artifacts aside
#   3. fetch data-index.tar.gz from the Release named in data-index.env into data/
#   4. python3 switch_source_path.py   -> rewrites `source:` to `path: ..`
#   5. conda-build conda.recipe
#   6. restore the URL recipe, delete data/index
nohup bash build_full_<VER>.sh > build_<VER>.log 2>&1 &
```
> `conda-build` lives in the **base** env (`/nfs/users/zhaojw/miniconda3/bin`).

Then:
```bash
export PATH=/nfs/users/zhaojw/miniconda3/bin:$PATH
sha256sum /nfs/users/zhaojw/miniconda3/conda-bld/noarch/mirdeep-p3-<VER>-0.conda
anaconda upload /nfs/users/zhaojw/miniconda3/conda-bld/noarch/mirdeep-p3-<VER>-0.conda
# no -u flag: `anaconda whoami` already reports "jaguares"
```
Verify:
```bash
anaconda show jaguares/mirdeep-p3                       # new version listed
conda create -n _verify -c jaguares -c conda-forge -c bioconda mirdeep-p3=<VER>
conda run -n _verify mirdeep-p3 --version               # MirDeep-P3 <VER>
```

## 7. bioconda

The recipe lives in `bioconda-recipes/recipes/mirdeep-p3` and is submitted from
the fork `mercedes-cykt/bioconda-recipes`, branch `add-mirdeep-p3`.

```bash
git clone --depth 1 --filter=blob:none --sparse \
  --branch add-mirdeep-p3 https://github.com/mercedes-cykt/bioconda-recipes.git bc-recipes
cd bc-recipes && git sparse-checkout set recipes/mirdeep-p3
# edit meta.yaml (version/url/sha256) and build.sh, then
git push origin add-mirdeep-p3     # updates the open PR
```

Non-obvious requirements for this recipe:

- **`data/index` must be a second source.**  The code resolves the databases as
  `<project_root>/data/index/*` with no CLI flag or environment variable to
  relocate them, so a package without them installs but cannot run
  `identification`.
  ```yaml
  source:
    - url: .../archive/refs/tags/<REL>.tar.gz
      sha256: <source sha256>
    - url: .../releases/download/mirdeep-p3-v3.1.4c-full/data-index.tar.gz
      sha256: <index sha256>
      folder: data-index
  ```
- **`folder:` still strips a single top-level directory.**  The index tarball's
  leading `index/` is removed, so the files land directly in
  `$SRC_DIR/data-index/`.  Probe for both layouts in `build.sh` and **fail
  loudly** if the index is missing — otherwise CI stays green while shipping a
  broken package.
- `about.license_file` needs the file **inside the recipe directory** — copy the
  repo `LICENSE` in.
- `PREFIX` is exported into the test script (`conda_build/build.py` sets
  `env["PREFIX"] = metadata.config.test_prefix`), so assertions such as
  `test -s "$PREFIX/share/mirdeep-p3/data/index/rfam_index.1.ebwt"` are valid and
  are what actually guards the bundled data.
- `noarch: generic` builds once (no per-platform fan-out), so the extra ~189 MB
  source download is paid once.

## 8. `data-index.*` — the load-bearing legacy Release

`DATA_INDEX_RELEASE` in `data-index.env` currently points at
`mirdeep-p3-v3.1.4c-full`.  That tag **and its Release must never be deleted**:
it hosts `data-index.tar.gz`, consumed by

- `data-index.env` -> `ci.yml`, `docker-release.yml`, `test.sh`
- the README / CONTRIBUTING download commands
- `conda.recipe` and the bioconda recipe (second `source:`)

Nothing about it follows the software version — that is intentional.  To move to
a new index, publish a new Release (a `data-v*` tag is a good name), change the
single line in `data-index.env`, and update the two `wget` lines in
README.md / CONTRIBUTING.md.  Do not delete the old Release until nothing
references it any more.

## 9. Credentials

| Credential | Identity | Can do |
|---|---|---|
| `~/.git-credentials` | `mercedes-cykt` | classic token, full `repo` scope — push to forks, comment on / edit PRs in **other** repos, delete refs |
| `~/.github_token` | `mercedes-cykt` | fine-grained PAT, **our repo only** — reading is fine, but writing to `bioconda-recipes` returns `403 Resource not accessible by personal access token` |

Use the classic one for anything outside `YangXZ-lab/mirdeep-p3`:
```bash
CRED=$(sed -n 's|^https://[^:]*:\([^@]*\)@github.com|\1|p' ~/.git-credentials | head -1)
```
Querying Actions:
```bash
curl -s -H "Authorization: token $CRED" \
  "https://api.github.com/repos/YangXZ-lab/mirdeep-p3/actions/runs?per_page=3"
```
`bin/mirdeep-p3` is stored in git as mode `100644`; the executable bit is applied
by `chmod 755` in CI and in the conda `build.sh`, so do not "fix" the mode.

## 10. Known flakiness and workarounds

| Symptom | Cause | Workaround |
|---|---|---|
| `git push` / `git fetch`: "Failed to connect to github.com port 443 after 13x s" | cluster egress to `github.com` is intermittent — TCP connect itself fails while `api.github.com` and `codeload.github.com` keep working | retry; use `api.github.com` for API work, `codeload` for downloads, or a pre-seeded cache |
| `conda-build`: `RuntimeError: Could not download <url>` | the same flakiness hitting conda's downloader (60 s read / 9 s connect, 3 retries) | pre-seed `<croot>/src_cache/` — see below |
| `docker pull` from compute1 is unusably slow | direct Docker Hub is blocked; only domestic mirrors work | use the GitHub Release image tarball (§5) |
| R: `ggsave(device = cairo_pdf)` writes nothing, no error | that R build | use base `pdf()` + `print()` |

Pre-seeding the conda-build source cache: the filename is
`<stem>_<sha256[:10]><ext>` — from `conda_build/source.py`:
`ext_re = re.compile(r"(.*?)(\.(?:tar\.)?[^.]+)$")` — stored in
`<croot>/src_cache/`.  With the right name in place the log says
`Found source in cache: ...` and nothing is downloaded.
```bash
C=/nfs/users/zhaojw/Genburten/temp/bc-bld/src_cache
cp /path/to/<REL>.tar.gz      $C/<REL>_<sha256[:10]>.tar.gz
cp /path/to/data-index.tar.gz $C/data-index_<sha256[:10]>.tar.gz
```

## 11. Post-release verification

- [ ] `docker run --rm merc3dez/mirdeep-p3:<VER>-full --version` -> `MirDeep-P3 <VER>`
- [ ] image contains `/opt/mirdeep-p3/data/index` (20 files, ~592 MB)
- [ ] `docker run --rm merc3dez/mirdeep-p3:<VER>-full -h` renders
- [ ] `conda create -c jaguares -c conda-forge -c bioconda mirdeep-p3=<VER>` then `mirdeep-p3 --version`
- [ ] Aliyun: `docker pull crpi-...cn-beijing.personal.cr.aliyuncs.com/merc3dez/mirdeep-p3:<VER>-full`
- [ ] bioconda PR CI green (`Lint`, `Linux Tests`, `OSX-64 Tests`, `build and test (ARM)`, `Summary`)
- [ ] GitHub Release `<REL>-full` has the image tarball attached; tag `<REL>` exists
- [ ] `data-index.tar.gz` still downloadable from `mirdeep-p3-v3.1.4c-full`
