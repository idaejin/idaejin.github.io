module github.com/idaejin/idaejin.github.io

go 1.19

require (
        github.com/wowchemy/wowchemy-hugo-themes/modules/wowchemy/v5 v5.0.0-00010101000000-000000000000
        github.com/wowchemy/wowchemy-hugo-themes/modules/wowchemy-plugin-netlify v1.0.0-00010101000000-000000000000
)

replace github.com/wowchemy/wowchemy-hugo-themes/modules/wowchemy/v5 => ./themes/github.com/wowchemy/wowchemy-hugo-themes/modules/wowchemy/v5

replace github.com/wowchemy/wowchemy-hugo-themes/modules/wowchemy-plugin-netlify => ./themes/github.com/wowchemy/wowchemy-hugo-themes/modules/wowchemy-plugin-netlify
