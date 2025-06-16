/^QUERY/ { print }
/Feature/ { print; n=2; next }
n > 0 { print; n-- }