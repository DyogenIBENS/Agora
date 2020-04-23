#!/bin/bash
for i in *.mmd
do
    [ "$i" -nt "${i/mmd/jpg}" ] && curl "https://mermaid.ink/img/$(base64 < "$i" | tr -d '\n')" > "${i/mmd/jpg}"
    #[ "$i" -nt "${i/mmd/svg}" ] && curl "https://mermaid.ink/svg/$(base64 < "$i" | tr -d '\n')" > "${i/mmd/svg}"
done
