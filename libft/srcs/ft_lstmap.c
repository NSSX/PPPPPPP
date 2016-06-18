/* ************************************************************************** */
/*                                                                            */
/*                                                        :::      ::::::::   */
/*   ft_lstmap.c                                        :+:      :+:    :+:   */
/*                                                    +:+ +:+         +:+     */
/*   By: igomez <igomez@student.42.fr>              +#+  +:+       +#+        */
/*                                                +#+#+#+#+#+   +#+           */
/*   Created: 2014/11/07 20:51:06 by igomez            #+#    #+#             */
/*   Updated: 2014/11/20 12:22:02 by igomez           ###   ########.fr       */
/*                                                                            */
/* ************************************************************************** */

#include "libft.h"

t_list			*ft_lstmap(t_list *lst, t_list *(*f)(t_list *elem))
{
	t_list	*res;
	t_list	*node;

	if (!lst || !f)
		return (NULL);
	res = NULL;
	while (lst)
	{
		if (!(node = ft_lstnew(lst->content, lst->content_size)))
			return (NULL);
		ft_lstadd_back(&res, f(node));
		lst = lst->next;
	}
	return (res);
}
