dnl Derive package version from git at autoconf time.
dnl Rules:
dnl - exact tag on HEAD -> that tag
dnl - otherwise, latest tag + "-dirty"
dnl - if no tags in repo -> "v0.1-dirty"
dnl - if not in git repo / git missing -> "UNKNOWN"
AC_DEFUN([CA_GIT_DESCRIBE_VERSION],
  [m4_esyscmd_s([sh -c '
if command -v git >/dev/null 2>&1 && git rev-parse --is-inside-work-tree >/dev/null 2>&1; then
  tag=$(git describe --tags --exact-match HEAD 2>/dev/null || true)
  if test -n "$tag"; then
    printf "%s" "$tag"
  else
    latest=$(git describe --tags --abbrev=0 2>/dev/null || true)
    if test -n "$latest"; then
      printf "%s" "${latest}-dirty"
    else
      printf "%s" "v0.1-dirty"
    fi
  fi
else
  printf "%s" "UNKNOWN"
fi
'])])
